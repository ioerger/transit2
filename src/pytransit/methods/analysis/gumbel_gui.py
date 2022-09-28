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

from pytransit.basics.lazy_dict import LazyDict

from pytransit.universal_data import universal
from pytransit.components.parameter_panel import panel as parameter_panel
from pytransit.components.parameter_panel import progress_update
from pytransit.components.panel_helpers import *
from pytransit.components.spreadsheet import SpreadSheet
import pytransit.tools.console_tools as console_tools
import pytransit.tools.gui_tools as gui_tools
import pytransit.tools.transit_tools as transit_tools
import pytransit.tools.tnseq_tools as tnseq_tools
import pytransit.tools.norm_tools as norm_tools
import pytransit.tools.stat_tools as stat_tools
import pytransit.basics.csv as csv
import pytransit.components.results_area as results_area
from pytransit.components.panel_helpers import make_panel, create_run_button, create_normalization_input, create_reference_condition_input, create_include_condition_list_input, create_exclude_condition_list_input, create_n_terminus_input, create_c_terminus_input, create_pseudocount_input, create_winsorize_input, create_alpha_input, create_button, create_text_box_getter, create_button, create_check_box_getter, create_control_condition_input, create_experimental_condition_input, create_preview_loess_button


command_name = sys.argv[0]

class Analysis:
    identifier  = "Gumbel"
    short_name  = "gumbel_gui"
    long_name   = "gumbel_gui"
    short_desc = "Bayesian analysis of essentiality based on long gaps."
    long_desc = """Bayesian methods of analyzing longest runs of non-insertions in a row. Estimates the parameters using the MCMC sampling, and estimates posterior probabilities of essentiality. 

    Reference: DeJesus et al. (2013; Bioinformatics)"""
    transposons = ["himar1"]
    columns = ["Orf", "Name", "Desc", "k", "n", "r", "s", "zbar", "Call"]
    
    inputs = LazyDict(
        combined_wig = None,
        metadata = None,
        condition = None, # all reps will be combined; later, allow user to select individual wigs files
        wig_files = None,
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
    # -n for normalization?
    ]

    usage_string = """python3 %s gumbel_gui <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -b <integer>    :=  Number of Burn-in samples. Default -b 500
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -t <integer>    :=  Trims all but every t-th value. Default: -t 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        -iN <float>     :=  Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 0
        """ % sys.argv[0]
    
    wxobj = None
    panel = None
    
    def __init__(self, *args, **kwargs):
        Analysis.instance = self
        self.full_name        = f"[{self.short_name}]  -  {self.short_desc}"
        self.transposons_text = transit_tools.get_transposons_text(self.transposons)
        self.filetypes        = [File]
        self.method           = Analysis # backwards compat
        self.gui              = self     # backwards compat
    
    def __str__(self):
        return f"""
            Analysis Method:
                Short Name:  {self.short_name}
                Long Name:   {self.long_name}
                Short Desc:  {self.short_desc}
                Long Desc:   {self.long_desc}
                GUI:         {self.gui}
        """.replace('\n            ','\n').strip()
    
    def __repr__(self):
        return f"{self.inputs}"

    def create_input_field(self, panel, sizer, label, value,tooltip=None):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text=label,
            default_value=value,
            tooltip_text=tooltip,
        )
        return lambda *args: get_text()

    def define_panel(self, _):
        from pytransit.tools.transit_tools import wx
        self.panel = make_panel()
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.value_getters = LazyDict()
        sample_getter                = create_text_box_getter(self.panel, main_sizer, label_text="Samples", default_value="10000", tooltip_text="")
        
        self.value_getters.condition       = create_condition_choice(self.panel,main_sizer,"Condition to analyze:")
        self.value_getters.normalization   = create_normalization_input(self.panel, main_sizer) # TTR 
        self.value_getters.samples         = lambda *args: int(sample_getter(*args))
        self.value_getters.burnin          = self.create_input_field(self.panel,main_sizer, label="Burnin",value=500,tooltip="Burnin")
        self.value_getters.trim            = self.create_input_field(self.panel,main_sizer, label="trim",value=1,tooltip="trim")
        self.value_getters.n_terminus      = create_n_terminus_input(self.panel, main_sizer)
        self.value_getters.c_terminus      = create_c_terminus_input(self.panel, main_sizer)


        create_run_button(self.panel, main_sizer)

        parameter_panel.set_panel(self.panel)
        self.panel.SetSizer(main_sizer)
        self.panel.Layout()
        main_sizer.Fit(self.panel)

    @classmethod
    def from_gui(cls, frame):
        with gui_tools.nice_error_log:
            # 
            # get wig files
            # 
            combined_wig = universal.session_data.combined_wigs[0]
            Analysis.inputs.combined_wig = combined_wig.main_path
            # assume all samples are in the same metadata file
            Analysis.inputs.metadata_path = universal.session_data.combined_wigs[0].metadata_path 


            
            # 
            # get annotation
            # 
            Analysis.inputs.annotation_path = universal.session_data.annotation_path
            # FIXME: enable this once I get a valid annotation file example
            # if not transit_tools.validate_annotation(Analysis.inputs.annotation):
            #     return None
            


            for each_key, each_getter in Analysis.instance.value_getters.items():
                try:
                    Analysis.inputs[each_key] = each_getter()
                except Exception as error:
                    raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
            ###transit_tools.log("included_conditions", Analysis.inputs.included_conditions)

            # 
            # get filename for results
            #Analysis.inputs.output_path = "%s.dat" % (Analysis.inputs.output_basename)

            Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
                default_file_name="output.dat",
                output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
            )
            # if not Analysis.inputs.output_path:
            #     return None

            # 
            # extract universal data
            # 
            # cwig_path     = universal.session_data.combined_wigs[0].main_path
            # metadata_path = universal.session_data.combined_wigs[0].metadata.path
            
            # Analysis.inputs.combined_wig_params = dict(
            #     combined_wig=cwig_path,
            #     samples_metadata=metadata_path,
            # )
            
            # Analysis.inputs.update(dict(
            #     annotation_path_exp=Analysis.inputs.annotation_path_exp
            # ))

            #if not Analysis.inputs.output_path: return None ### why?
            return Analysis.instance

    @classmethod
    def from_args(cls, args, kwargs): # clean_args() was already called in pytransit/__main__.py

      if len(args)!=3: print(cls.usage_string); sys.exit(0) # use transit_error()?

      # check for unrecognized flags
      console_tools.handle_unrecognized_flags(
            "-s -b -m -t -r -iN -iC".split(),
            kwargs,
            cls.usage_string,
      )

      Analysis.inputs.update(dict(
        combined_wig = None,
        metadata = None,
        wig_files = args[0].split(','),
        annotation_path = args[1],
        output_path = args[2],
        normalization = kwargs.get("n", "TTR"),
        samples = int(kwargs.get("s", 10000)),
        burnin = int(kwargs.get("b", 500)),
        read_count = int(kwargs.get("r", 1)),
        trim = int(kwargs.get("t", 1)),
        replicates = kwargs.get("r", "Sum"),
        iN = float(kwargs.get("iN", 0.00)), 
        iC=float(kwargs.get("iC", 0.00)),
      ))

       
        
      return Analysis.instance
        
    def Run(self):
 
        with gui_tools.nice_error_log:
            transit_tools.log("Starting gumbel analysis")
            self.start_time = time.time()

            #######################
            # get data

            if self.inputs.combined_wig!=None:  # assume metadata and condition are defined too
              transit_tools.log("Getting Data from %s" % self.inputs.combined_wig)
              position, data, filenames_in_comb_wig = tnseq_tools.read_combined_wig(self.inputs.combined_wig)

              metadata = tnseq_tools.CombinedWigMetadata(self.inputs.metadata_path)
              indexes = {}
              for i,row in enumerate(metadata.rows): 
                cond = row["Condition"] 
                if cond not in indexes: indexes[cond] = []
                indexes[cond].append(i)
              cond = Analysis.inputs.condition
              ids = [metadata.rows[i]["Id"] for i in indexes[cond]]
              transit_tools.log("selected samples for gumbel (cond=%s): %s" % (cond,','.join(ids)))
              data = data[indexes[cond]] # project array down to samples selected by condition

              # now, select the columns in data corresponding to samples that are replicates of desired condition...

            elif self.inputs.wig_files!=None:
              transit_tools.log("Getting Data")
              (data, position) = transit_tools.get_validated_data( self.inputs.wig_files, wxobj=self.wxobj )

            else: print("error: must provide either combined_wig or list of wig files"); sys.exit(0) ##### use transit.error()?

            (K, N) = data.shape
            merged = numpy.sum(data, axis=0)
            self.inputs.nsites, nzeros = (
                merged.shape[0],
                numpy.sum(merged == 0),
            )  # perhaps I should say >minCount
            self.inputs.sat = (self.inputs.nsites - nzeros) / float(self.inputs.nsites)

            # normalize the counts
            if self.inputs.normalization and self.inputs.normalization != "nonorm":
                transit_tools.log("Normalizing using: %s" % self.inputs.normalization)
                (data, factors) = norm_tools.normalize_data( data, self.inputs.normalization, self.inputs.wig_files, self.inputs.annotation_path )
                
                # read-in genes from annotation
                G = tnseq_tools.Genes(
                    self.inputs.wig_files,
                    self.inputs.annotation_path,
                    data=data,
                    position=position,
                    minread=1,  ### add these options?
                    reps=self.inputs.replicates,
                    #ignore_codon=self.ignore_codon,
                    n_terminus=self.inputs.iN, 
                    c_terminus=self.inputs.iC,
                )
                N = len(G)

            # could also read-in gumbel_results_file as csv here

            ###########################
            # process data

            transit_tools.log("processing data")
            Z_sample, phi_sample, count, acctot= self.calc_gumbel(G)

            ###########################
            # write output        

            self.write_gumbel_results(G, Z_sample, phi_sample, count, acctot)
            results_area.add(self.inputs.output_path)
            transit_tools.log(f"Finished running {Analysis.short_name}")       

    def calc_gumbel(self,G):
        transit_tools.log("Starting Gumbel Method")

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
        transit_tools.log("Reading Annotation")
        ii_good = numpy.array([self.good_orf(g) for g in G])  # Gets index of the genes that can be analyzed

        K = G.local_insertions()[ii_good]
        N = G.local_sites()[ii_good]
        R = G.local_runs()[ii_good]
        S = G.local_gap_span()[ii_good]
        T = G.local_gene_span()[ii_good]

        ####################################################
        transit_tools.log("Doing Regression")
        mu_s, temp, sigma_s = stat_tools.regress(R, S)  # Linear regression to estimate mu_s, sigma_s for span data
        mu_r, temp, sigma_r = stat_tools.regress(S, R)  # Linear regression to estimate mu_r, sigma_r for run data

        N_GENES = len(G)
        N_GOOD = sum(ii_good)

        transit_tools.log("Setting Initial Class")
        Z_sample = numpy.zeros((N_GOOD, self.samples))
        Z = [self.classify(g.n, g.r, 0.5) for g in G if self.good_orf(g)]
        Z_sample[:, 0] = Z
        N_ESS = numpy.sum(Z_sample[:, 0] == 1)

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
        while i < self.samples:

            try:
                # PHI
                acc = 1.0
                phi_new = phi_old + random.gauss(mu_c, sigma_c)
                i0 = Z_sample[:, i - 1] == 0
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
                    Z_sample[:, i] = Z
                    i += 1

            except ValueError as e:
                transit_tools.log("Error: %s" % e)
                transit_tools.log("This is likely to have been caused by poor data (e.g. too sparse)." )
                transit_tools.log("If the density of the dataset is too low, the Gumbel method will not work.")
                transit_tools.log("Quitting.")
                return

            #            print(i,phi_new,w1,G[idxG].name,N[idxN],R[idxN],Z[idxN])

            phi_old = phi_new
            # Update progress
            percentage = (100.0 * (count + 1) / (self.samples + self.burnin))
            text = "Running Gumbel Method with Binomial Essentiality Calls... %5.1f%%" % percentage
            progress_update(text, percentage)


        return (Z_sample, phi_sample, count, acctot)

    def write_gumbel_results(self, G, Z_sample, phi_sample, count, acctot):
        
        ZBAR = numpy.apply_along_axis(numpy.mean, 1, Z_sample)
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

            # data.append(
            #     "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%s\n"
            #     % (g.orf, g.name, g.desc, g.k, g.n, g.r, g.s, zbar, call)
            # )
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
            file_kind=Analysis.identifier,
            rows=rows,
            column_names=Analysis.columns,
            extra_info=dict(
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
                time=(time.time() - self.start_time),

                ES = str(calls.count("E")) + " #essential based on Gumbel",
                ESB = str(calls.count("EB")) + " #essential based on Binomial",
                NE = str(calls.count("NE")) + " #non-essential",
                U = str(calls.count("U")) + " #uncertain",
                S = str(calls.count("S")) +" #too-short",
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
class File(Analysis):
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Analysis.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=transit_tools.basename(self.path),
            type=Analysis.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title=Analysis.identifier,heading="",column_names=self.column_names,rows=self.rows).Show(),
            })
        )
        
        # 
        # get column names
        # 
        comments, headers, rows = csv.read(self.path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
        if len(comments) == 0:
            raise Exception(f'''No comments in file, and I expected the last comment to be the column names, while to load tnseq_stats file "{self.path}"''')
        self.column_names = comments[-1].split("\t")
        
        # 
        # get rows
        #
        self.rows = []
        for each_row in rows:
            row = {}
            for each_column_name, each_cell in zip(self.column_names, each_row):
               row[each_column_name] = each_cell
            self.rows.append(row)
        
    
    def __str__(self):
        return f"""
            File for {self.short_name}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
    
    
Method = GUI = Analysis
Analysis() # make sure there's one instance
