import sys
import os
import time
import ntpath
import math
import random
import datetime
import itertools
import statistics
import heapq
import collections
import numpy
import scipy

from pytransit.generic_tools.lazy_dict import LazyDict

from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.components.parameter_panel import progress_update, set_instructions
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.generic_tools import csv, misc, informative_iterator
from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools, stat_tools
import pytransit.components.results_area as results_area

@misc.singleton
class Method:
    name = "Gumbel"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name} - Bayesian analysis of essentiality based on long gaps."
    description = f"""Bayesian methods of analyzing longest runs of non-insertions in a row. Estimates the parameters using the MCMC sampling, and estimates posterior probabilities of essentiality. 

    Reference: DeJesus et al. (2013; Bioinformatics)"""
    
    transposons = ["himar1"]
    column_names = [
        "ORF", 
        "Gene Name", 
        "Description", 
        "Number Of Insertions Within ORF",
        "Total Number Of TA Sites Within ORF", 
        "Length Of Maximum Run Of Non Insertions",
        "Nucleotide Span For Maximum Run Of Non Insertions",
        "Posterior Probability Of Essentiality", # Z Bar
        "Essentiality Call",
    ]
    
    inputs = LazyDict(
        combined_wig = None,
        metadata = None,
        condition = None, # all reps will be combined; later, allow user to select individual wigs files
        annotation_path = None,
        output_path = None,
        normalization = "TTR",
        samples = 10000,
        burnin = 500,
        read_count = 1,
        trim = 1,
        replicates = "Sum",
        iN = 0,
        iC=0,

        cache_expruns = {},
        cache_nn = {},

        EXACT = 20,
        ALPHA = 1,
        BETA =1,
    )
    
    valid_cli_flags = [
        "--s",
        "--b",
        "--m",
        "--t",
        "--r",
        "--iN",
        "--iC",
    ]

    usage_string = f"""
        Usage:
            {console_tools.subcommand_prefix} gumbel <combined_wig_file> <metadata_file> <annotation_file> <condition_to_analyze> <output_file> [Optional Arguments]
    
        Optional Arguments:
            --s <integer> := Number of samples. Default: --s 10000
            --b <integer> := Number of Burn-in samples. Default --b 500
            --m <integer> := Smallest read-count to consider. Default: --m 1
            --t <integer> := Trims all but every t-th value. Default: --t 1
            --r <string>  := How to handle replicates. Sum or Mean. Default: --r Sum
            --iN <float>  := Ignore TAs occurring within given percentage (as integer) of the N terminus. Default: --iN 0
            --iC <float>  := Ignore TAs occurring within given percentage (as integer) of the C terminus. Default: --iC 0
    """.replace("\n        ", "\n")
    
    @gui.add_menu("Method", "Himar1", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            set_instructions(
                title_text=self.name,
                sub_text="",
                method_specific_instructions="""
                    The Gumbel can be used to determine which genes are essential in a single condition. It does a gene-by-gene analysis of the insertions at TA sites with each gene, makes a call based on the longest consecutive sequence of TA sites without insertion in the genes, calculates  the probability of this using a Bayesian model.
                        
                    1. Of the Conditions in the Conditions pane, select one using the 'Conditions' dropdown

                    2. [Optional] Select/Adjust values for the remaining parameters

                    3. Click Run
                """.replace("\n                    ","\n"),
            )
                
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            self.value_getters = LazyDict()
            
            self.value_getters.condition       = panel_helpers.create_condition_input(panel, main_sizer)
            self.value_getters.normalization   = panel_helpers.create_normalization_input(panel, main_sizer) # TTR 
            self.value_getters.samples         = panel_helpers.create_int_getter(panel, main_sizer, label_text="Samples", default_value=10000, tooltip_text="The default setting is to run the simulation for 10,000 iterations. This is usually enough to assure convergence of the sampler and to provide accurate estimates of posterior probabilities. Less iterations may work, but at the risk of lower accuracy.")
            self.value_getters.burnin          = panel_helpers.create_int_getter(panel, main_sizer, label_text="Burnin", default_value=500, tooltip_text="Because the MH sampler many not have stabilized in the first few iterations, a “burn-in” period is defined. Samples obtained in this “burn-in” period are discarded, and do not count towards estimates.")
            self.value_getters.trim            = panel_helpers.create_int_getter(panel, main_sizer, label_text="Trim", default_value=1, tooltip_text="The MH sampler produces Markov samples that are correlated. This parameter dictates how many samples must be attempted for every sampled obtained. Increasing this parameter will decrease the auto-correlation, at the cost of dramatically increasing the run-time. For most situations, this parameter should be left at the default of “1”.")
            self.value_getters.n_terminus      = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus      = panel_helpers.create_c_terminus_input(panel, main_sizer)
            self.value_getters.read_count      = panel_helpers.create_int_getter(panel, main_sizer, label_text="Minimum Read Count", default_value=1, tooltip_text="The minimum read count that is considered a true read. Because the Gumbel method depends on determining gaps of TA sites lacking insertions, it may be susceptible to spurious reads (e.g. errors). The default value of 1 will consider all reads as true reads. A value of 2, for example, will ignore read counts of 1.")
            self.value_getters.replicates      = panel_helpers.create_choice_input(panel, main_sizer, label="Replicates:", options=["Sum", "Mean", "TTRMean"], tooltip_text="Determines how to handle replicates, and their read-counts. When using many replicates, using 'Mean' may be recommended over 'Sum'")
        
    
    @staticmethod
    def from_gui(frame):
        with gui_tools.nice_error_log:
            # 
            # get wig files
            # 
            combined_wig = gui.combined_wigs[-1]
            Method.inputs.combined_wig_path = combined_wig.main_path
            # assume all samples are in the same metadata file
            Method.inputs.metadata_path = gui.combined_wigs[-1].metadata_path 


            
            # 
            # get annotation
            # 
            Method.inputs.annotation_path = gui.annotation_path

            print(Method.inputs)

            for each_key, each_getter in Method.value_getters.items():
                try:
                    Method.inputs[each_key] = each_getter()
                except Exception as error:
                    raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')

            Method.inputs.output_path = gui_tools.ask_for_output_file_path(
                default_file_name=f"{Method.cli_name}_output.tsv",
                output_extensions=transit_tools.result_output_extensions,
            )

            # 
            # validate
            # 
            assert Method.inputs.condition != "[None]", "Please select a condition"

            #if not Method.inputs.output_path: return None ### why?
            return Method

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=5)

        Method.inputs.update(dict(
            combined_wig_path=args[0],
            metadata_path=args[1],
            annotation_path=args[2],
            condition=args[3],
            output_path=args[4],
            normalization=kwargs.get("n", "TTR"),
            samples=int(kwargs.get("s", 10000)),
            burnin=int(kwargs.get("b", 500)),
            read_count=int(kwargs.get("m", 1)),
            trim=int(kwargs.get("t", 1)),
            replicates=kwargs.get("r", "Sum"),
            iN=float(kwargs.get("iN", 0.00)), 
            iC=float(kwargs.get("iC", 0.00)),
        ))

        Method.Run()
        
    def Run(self):
        with gui_tools.nice_error_log:
            logging.log("Starting gumbel analysis")
            self.start_time = time.time()
            logging.log("Processing data from %s" % self.inputs.combined_wig_path)
            
            combined_wig = tnseq_tools.CombinedWig(
                main_path=self.inputs.combined_wig_path,
                metadata_path=self.inputs.metadata_path,
                annotation_path=self.inputs.annotation_path,
            )
            
            # make sure selected condition exists
            if self.inputs.condition not in combined_wig.condition_names:
                logging.error(f"\n\nThe condition name {self.inputs.condition} is not one of the available conditions: {combined_wig.condition_names}\n")
            
            # create combined wig with just the selected condition, normalize it, then organize the data by-genes
            normalized_combined_wig = tnseq_tools.CombinedWig(
                main_path=self.inputs.combined_wig_path,
                metadata_path=self.inputs.metadata_path,
                annotation_path=self.inputs.annotation_path,
            ).with_only(
                condition_names=[self.inputs.condition]
            ).normalized_with( # good to normalize AFTER filtering out the other conditions
                kind=self.inputs.normalization
            )
            
            genes_object = normalized_combined_wig.get_genes(
                n_terminus=self.inputs.iN,
                c_terminus=self.inputs.iC,
                reps=self.inputs.replicates,
                minread=1,   ### add these options?
            )
            
            logging.log("Computing Gumbel")
            z_sample, phi_sample, count, acctot = self.calc_gumbel(genes_object)
            self.write_gumbel_results(genes_object, z_sample, phi_sample, count, acctot)
            results_area.add(self.inputs.output_path)
            logging.log(f"Finished running {Method.identifier}")       

    def calc_gumbel(self,G):
        logging.log("Starting Gumbel Method")

        # Set Default parameter values
        w1 = 0.15
        w0 = 1.0 - w1
        self.ALPHA = 1
        self.BETA =1
        self.ALPHA_w = 600
        BETA_w = 3400
        mu_c = 0
        acctot = 0.0
        phi_start = 0.3
        sigma_c = 0.01
        self.samples = 10000
        self.burnin = 500
        self.trim=1

        # Get orf data
        logging.log("Reading Annotation")
        ii_good = numpy.array([self.good_orf(g) for g in G])  # Gets index of the genes that can be analyzed

        K = G.local_insertions()[ii_good]
        N = G.local_sites()[ii_good]
        R = G.local_runs()[ii_good]
        S = G.local_gap_span()[ii_good]
        T = G.local_gene_span()[ii_good]

        ####################################################
        logging.log("Doing Regression")
        mu_s, temp, sigma_s = stat_tools.regress(R, S)  # Linear regression to estimate mu_s, sigma_s for span data
        mu_r, temp, sigma_r = stat_tools.regress(S, R)  # Linear regression to estimate mu_r, sigma_r for run data

        N_GENES = len(G)
        N_GOOD = sum(ii_good)

        logging.log("Setting Initial Class")
        z_sample = numpy.zeros((N_GOOD, self.samples))
        Z = [self.classify(g.n, g.r, 0.5) for g in G if self.good_orf(g)]
        z_sample[:, 0] = Z
        N_ESS = numpy.sum(z_sample[:, 0] == 1)

        phi_sample = numpy.zeros(self.samples)  # []
        phi_sample[0] = phi_start
        phi_old = phi_start
        phi_new = 0.00

        SIG = numpy.array(
            [
                self.sigmoid(g.s, g.t) * scipy.stats.norm.pdf(g.r, mu_r * g.s, sigma_r)
                for g in G
                if self.good_orf(g)
            ]
        )

        i = 1
        count = 0
        cli_progress_bar = iter(
            informative_iterator.ProgressBar(
                self.samples + self.burnin,
                title="Running Gumbel ..."
            )
        )
        while i < self.samples:
            progress, _ = next(cli_progress_bar)
            try:
                # PHI
                acc = 1.0
                phi_new = phi_old + random.gauss(mu_c, sigma_c)
                i0 = z_sample[:, i - 1] == 0
                if (
                    phi_new > 1
                    or phi_new <= 0
                    or (
                        self.f_non(phi_new, N[i0], R[i0])
                        - self.f_non(phi_old, N[i0], R[i0])
                    )
                    < math.log(random.uniform(0, 1))
                ):
                    phi_new = phi_old
                    acc = 0.0
                    flag = 0

                # Z
                Z = self.sample_Z(phi_new, w1, N, R, S, T, mu_s, sigma_s, SIG)

                # w1
                N_ESS = sum(Z == 1)
                w1 = scipy.stats.beta.rvs(N_ESS + self.ALPHA_w, N_GOOD - N_ESS + BETA_w)

                count += 1
                acctot += acc

                if (count > self.burnin) and (count % self.trim == 0):
                    phi_sample[i] = phi_new
                    z_sample[:, i] = Z
                    i += 1

            except ValueError as error:
                logging.log("Error: %s" % error)
                logging.log("This is likely to have been caused by poor data (e.g. too sparse)." )
                logging.log("If the density of the dataset is too low, the Gumbel method will not work.")
                logging.log("Quitting.")
                return

            phi_old = phi_new
            # Update progress
            percentage = (100.0 * (count + 1) / (self.samples + self.burnin))
            if gui.is_active:
                text = "Running Gumbel... %5.1f%%" % percentage
                progress_update(text, percentage)

        print() # to clear out iterator that may not finish (because of while-loop instead of for-loop)
        
        return (z_sample, phi_sample, count, acctot)

    def write_gumbel_results(self, G, z_sample, phi_sample, count, acctot):
        
        ZBAR = numpy.apply_along_axis(numpy.mean, 1, z_sample)
        (ess_t, non_t) = stat_tools.bayesian_essentiality_thresholds(ZBAR)
        binomial_n = math.log10(0.05) / math.log10(G.global_phi())

        i = 0
        data, calls = [], []
        for j, g in enumerate(G):
            if not self.good_orf(g):
                zbar = -1.0
            else:
                zbar = ZBAR[i]
                i += 1
            if zbar > ess_t:
                call = "E"
            elif G.local_sites()[j] > binomial_n and G.local_thetas()[j] == 0.0:
                call = "EB"
            elif non_t <= zbar <= ess_t:
                call = "U"
            elif 0 <= zbar < non_t:
                call = "NE"
            else:
                call = "S"
            data.append([g.orf, g.name, g.desc, g.k, g.n, g.r, g.s, zbar, call])
            calls.append(call)
        data.sort()


        rows = []
        for row_index, row in enumerate(data):
            (orf,name,desc, k,n, r, s, zbar,Call,) = row
            rows.append(("%s\t%s\t%s\t%d\t%d\t%d\t%1.2f\t%1.1f\t%s" % (orf,name,desc, k,n, r, s, zbar,Call)).split('\t'))

        # 
        # write to file
        # 
        transit_tools.write_result(
            path=self.inputs.output_path,
            file_kind=Method.identifier,
            rows=rows,
            column_names=Method.column_names,
            extra_info=dict(
                calculation_time=f"{(time.time() - self.start_time):0.1f}seconds",
                analysis_type=Method.identifier,
                files=dict(
                    combined_wig=Method.inputs.combined_wig_path,
                    annotation_path=Method.inputs.annotation_path,
                ),
                parameters=dict(
                    samples=self.inputs.samples,
                    norm=self.inputs.normalization,
                    burnin = self.inputs.burnin,
                    read_count = self.inputs.read_count,
                    trim = self.inputs.trim,
                    replicates = self.inputs.replicates,
                    iN = self.inputs.iN,
                    iC=self.inputs.iC,
                ),
                annotation_path=self.inputs.annotation_path,
                summary_info=dict(
                    ES = str(calls.count("E")),
                    ESB = str(calls.count("EB")),
                    NE = str(calls.count("NE")),
                    U = str(calls.count("U")),
                    S = str(calls.count("S")),
                ),
                naming_reference = dict(
                    ES = "Essential based on Gumbel",
                    ESB = "Essential based on Binomial",
                    NE = "Non-essential",
                    U = "Uncertain",
                    S = "Too-short",
                ),
            ),
        )
        
        
         

    def good_orf(self, gene):
        return gene.n >= 3 and gene.t >= 150

    def expected_runs_cached(self, n, q):
        if (n, q) not in self.inputs.cache_expruns:
            self.inputs.cache_expruns[(n, q)] = tnseq_tools.expected_runs(n, q)
        return self.inputs.cache_expruns[(n, q)]

    def classify(self, n, r, p):
        if n == 0:
            return 0
        q = 1 - p
        B = 1 / math.log(1 / p)
        u = math.log(n * q, 1 / p)
        BetaGamma = B * tnseq_tools.get_gamma()
        if (
            n < self.inputs.EXACT
        ):  # estimate more accurately based on expected run len, using self.EXACT calc for small genes
            exprun = self.expected_runs_cached(n, p)
            u = (
                exprun - BetaGamma
            )  # u is mu of Gumbel (mean=mu+gamma*beta); matching of moments
            # https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/
        pval = 1 - numpy.exp(scipy.stats.gumbel_r.logcdf(r, u, B))
        if pval < 0.05:
            return 1
        else:
            return 0

    def f_non(self, p, N, R):  # pass in P_nonins as p
        q = 1.0 - p
        BetaGamma = tnseq_tools.get_gamma() / math.log(1 / p)
        total = numpy.log(scipy.stats.beta.pdf(p, self.ALPHA, self.BETA))
        mu = numpy.log(N * q) / numpy.log(1 / p)
        for i in range(
            len(N)
        ):  # estimate more accurately based on expected run len, using self.EXACT calc for small genes
            if N[i] < self.inputs.EXACT:
                mu[i] = self.expected_runs_cached(int(N[i]), p) - BetaGamma
        sigma = 1 / math.log(1 / p)
        # for i in range(len(N)): print('\t'.join([str(x) for x in N[i],R[i],self.expected_runs_cached(int(N[i]),q),mu[i],scipy.stats.gumbel_r.pdf(R[i], mu[i], sigma)]))
        total += numpy.sum(scipy.stats.gumbel_r.logpdf(R, mu, sigma))
        return total

    def sample_Z(self, p, w1, N, R, S, T, mu_s, sigma_s, SIG):
        G = len(N)
        q = 1.0 - p
        BetaGamma = tnseq_tools.get_gamma() / math.log(1 / p)
        mu = numpy.log(N * q) / numpy.log(1 / p)
        for i in range(
            len(N)
        ):  # estimate more accurately based on expected run len, using self.EXACT calc for small genes
            if N[i] < self.inputs.EXACT:
                mu[i] = self.expected_runs_cached(int(N[i]), p) - BetaGamma
        sigma = 1.0 / math.log(1.0 / p)
        h0 = (
            (numpy.exp(scipy.stats.gumbel_r.logpdf(R, mu, sigma)))
            * scipy.stats.norm.pdf(S, mu_s * R, sigma_s)
            * (1 - w1)
        )
        h1 = SIG * w1
        h1 += 1e-10
        h0 += 1e-10  # to prevent div-by-zero; if neither class is probable, p(z1) should be ~0.5
        p_z1 = h1 / (h0 + h1)
        return scipy.stats.binom.rvs(1, p_z1, size=G)

    def sigmoid(self, d, n):
        Kn = 0.1
        MEAN_DOMAIN_SPAN = 300

        if d == 0:
            return 0.00
        f = 1.0 / (1.0 + math.exp(Kn * (MEAN_DOMAIN_SPAN - d)))
        # if n in self.cache_nn: return f/self.cache_nn[n]
        tot = 0
        N = int(n + 1)
        for i in range(1, N):
            tot += 1.0 / (1.0 + math.exp(Kn * (MEAN_DOMAIN_SPAN - i)))
        self.inputs.cache_nn[n] = tot
        return f / tot


@transit_tools.ResultsFile
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=transit_tools.basename(self.path),
            type=Method.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title=Method.identifier,heading=misc.human_readable_data(self.extra_data),column_names=self.column_names,rows=self.rows).Show(),
            })
        )
        
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        summary = self.extra_data.get("summary_info", {})
        summary_str = [str(summary[key])+" "+str(key) for key in ["ES","ESB", "NE","S","U"]] 
        self.values_for_result_table.update({"summary": "; ".join(summary_str) })
        
        parameters = self.extra_data.get("parameters",{})
        parameters_str = [str(key)+" : "+str(parameters[key]) for key in ["replicates","norm"]]
        self.values_for_result_table.update({"parameters": "; ".join(parameters_str) })


    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
    
    
