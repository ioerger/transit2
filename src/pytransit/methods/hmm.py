import sys
import os
import time
import ntpath
import math
import random
import datetime
import heapq
import collections

import numpy
import scipy
import scipy.stats
import heapq
import math
from pytransit.generic_tools.lazy_dict import LazyDict

from pytransit.specific_tools.transit_tools import wx, basename
from pytransit.specific_tools.tnseq_tools import Wig
from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.components.parameter_panel import progress_update, set_instructions
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.generic_tools import csv
import pytransit.generic_tools.misc as misc


    
@misc.singleton
class Method:
    name = "HMM"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name} - Method of genomic regions using a Hidden Markov Model"
    description = """Method of essentiality in the entire genome using a Hidden Markov Model. Capable of determining regions with different levels of essentiality representing Essential, Growth-Defect, Non-Essential and Growth-Advantage regions. Reference: DeJesus et al. (2013; BMC Bioinformatics)"""
    
    transposons = ["himar1"]
    categories = ["ES", "NE", "GD", "GA"]
    inputs = LazyDict(
        data_sources=[],
        ctrl_read_counts=None,
        ctrl_positions=None,
        
        annotation_path=None,
        output_path=None,
        replicates="Mean",
        normalization="TTR",
        loess_correction=False,
        n_terminus=0.0,
        c_terminus=0.0,
        conf_on=True,
    )
    
    valid_cli_flags = [
        "--r",
        "--n",
        "-l",
        "--iN",
        "--iC",
        "-conf-on",
    ]
    
    usage_string = f"""
        Usage:
            {console_tools.subcommand_prefix} hmm <combined_wig_file> <metadata_file> <annotation_file> <condition_to_analyze> <output_file>

        Optional Arguments:
            --r <string>     :=  How to handle replicates. Sum, Mean. Default: --r Mean
            --n <string>     :=  Normalization method. Default: --n TTR
            -l               :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Off.
            --iN <float>     :=  Ignore TAs occurring within given percentage (as integer) of the N terminus. Default: --iN 0
            --iC <float>     :=  Ignore TAs occurring within given percentage (as integer) of the C terminus. Default: --iC 0
            -conf-on         :=  enable additional columns with confidence information
    """.replace("\n        ", "\n")
    
    column_names = [
        "ORF",
        "Gene",
        "Annotation",
        "TAs",
        "ES Sites",
        "GD Sites",
        "NE Sites",
        "GA Sites",
        "Saturation",
        "Mean",
        "Call",
    ]
    
    @gui.add_menu("Method", "Himar1", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            set_instructions(
                title_text= self.name,
                sub_text= "Hidden Markov Model",
                method_specific_instructions="""
                    The HMM method can be used to determine the essentiality of the entire genome, as opposed to gene-level analysis of the other methods. It is capable of identifying regions that have unusually high or unusually low read counts (i.e. growth advantage or growth defect regions), in addition to the more common categories of essential and non-essential.
                    
                    1. Select a condition from the conditions panel

                    2. [Optional] Select/Adjust other parameters

                    3. Click Run
                """.replace("\n                    ","\n"),
            )
            # 
            # parameter inputs
            # 
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            self.value_getters = LazyDict()
            self.value_getters.condition              = panel_helpers.create_condition_input(panel, main_sizer)
            self.value_getters.normalization          = panel_helpers.create_normalization_input(panel, main_sizer)
            self.value_getters.replicates             = panel_helpers.create_choice_input(panel, main_sizer, label="Replicates:", options=["Mean", "Sum", "TTRMean"], tooltip_text="Determines how to handle replicates, and their read-counts. When using many replicates, using 'Mean' may be recommended over 'Sum'")
            self.value_getters.n_terminus             = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus             = panel_helpers.create_c_terminus_input(panel, main_sizer)
            self.value_getters.loess_correction       = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="LOESS Correction", default_value=False, tooltip_text="This adjusts insertion counts to correct for genome positional bias. Check to correct read-counts for possible regional biase using LOESS correction. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
            self.value_getters.conf_on                = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="Add confidence info columns", default_value=True, tooltip_text="")
        
    @staticmethod
    def from_gui(frame):
        # 
        # get wig files
        # 
        combined_wig = gui.combined_wigs[-1]
        Method.inputs.combined_wig = combined_wig.main_path
        Method.inputs.metadata     = combined_wig.metadata.path
        
        # 
        # get annotation
        # 
        Method.inputs.annotation_path = gui.annotation_path
        
        # 
        # setup custom inputs
        # 
        for each_key, each_getter in Method.value_getters.items():
            try:
                Method.inputs[each_key] = each_getter()
            except Exception as error:
                raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
        
        # 
        # validate
        # 
        assert Method.inputs.condition != "[None]", "Please select a condition"
        
        # 
        # save result files
        # 
        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}_output.tsv",
            output_extensions=transit_tools.result_output_extensions,
        )
        if not Method.inputs.output_path:
            return None
        
        Method.inputs.data_sources = [ gui.combined_wigs[-1].main_path ]
        
        # 
        # extract universal data
        # 
        Method.inputs.ctrl_read_counts, Method.inputs.ctrl_positions = transit_tools.gather_sample_data_for(conditions=[ Method.value_getters.condition() ])
        
        return Method

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=5)
        
        combined_wig    = args[0]
        metadata_path   = args[1]
        annotation_path = args[2]
        condition       = args[3]
        output_path     = args[4]
        
        replicates    = kwargs.get("r", Method.inputs.replicates)
        normalization = kwargs.get("n", Method.inputs.normalization)
        loess_correction         = kwargs.get("l", Method.inputs.loess_correction)
        n_terminus    = float(kwargs.get("iN", Method.inputs.n_terminus))
        c_terminus    = float(kwargs.get("iC", Method.inputs.c_terminus))
        conf_on      = float(kwargs.get("conf-on", False))

        ##################
        # read data      
        # old way: (ctrl_read_counts, ctrl_positions) = transit_tools.get_validated_data(ctrldata) 
        # read counts from combined_wig, like gumbel
        logging.log("Getting Data from %s" % combined_wig)
        position, data, filenames_in_comb_wig = tnseq_tools.CombinedWigData.load(combined_wig)
        metadata = tnseq_tools.CombinedWigMetadata(metadata_path)

        # now, select the columns in data corresponding to samples that are replicates of desired condition...
        indexes,ids = [],[]
        for i,f in enumerate(filenames_in_comb_wig): 
           cond = metadata.conditions_by_wig_fingerprint.get(f, "FLAG-UNMAPPED-CONDITION-IN-WIG")
           id = metadata.id_for(f)
           if cond==condition:
             indexes.append(i)
             ids.append(id)

        logging.log("selected samples for HMM analysis (cond=%s): %s" % (condition,','.join(ids)))
        data = data[indexes]
        if len(data)==0: print("error: it looks like there is a problem selecting samples"); sys.exit(0) # should I use error logging or transit_error?
        ctrldata = combined_wig
        ctrl_read_counts,ctrl_positions = data,position
        ##################

        Method.inputs.update(dict(
            combined_wig=combined_wig,
            metadata_path=metadata_path,
            annotation_path=annotation_path,
            condition=condition,
            output_path=output_path,
            data_sources=ctrldata,
            ctrl_read_counts=ctrl_read_counts,
            ctrl_positions=ctrl_positions,
            replicates=replicates,
            normalization=normalization,
            loess_correction=loess_correction,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            conf_on=conf_on,
        ))
        Method.Run()

    def Run(self):
        with gui_tools.nice_error_log:
            import pytransit.components.results_area as results_area
            # 
            # Calculations
            # 
            if True:
                self.max_iterations = len(Method.inputs.ctrl_positions) * len(Method.categories) + 1
                self.count = 1
                logging.log("Starting HMM Method")
                start_time = time.time()
                
                data, position = Method.inputs.ctrl_read_counts, Method.inputs.ctrl_positions
                (K, N) = data.shape

                # Normalize data
                if self.inputs.normalization != "nonorm":
                    logging.log("Normalizing using: %s" % self.inputs.normalization)
                    (data, factors) = norm_tools.normalize_data(
                        data, self.inputs.normalization, annotation_path=self.inputs.annotation_path,
                    )

                # Do loess_correction
                if self.inputs.loess_correction:
                    logging.log("Performing loess_correction Correction")
                    from pytransit.specific_tools import stat_tools
                    for j in range(K):
                        data[j] = stat_tools.loess_correction(position, data[j])

                hash = transit_tools.get_pos_hash(self.inputs.annotation_path)
                rv2info = tnseq_tools.AnnotationFile(path=self.inputs.annotation_path).orf_to_info

                if len(self.inputs.ctrl_read_counts) > 1:
                    logging.log("Combining Replicates as '%s'" % self.inputs.replicates)
                O = tnseq_tools.combine_replicates(data, method=self.inputs.replicates) + 1
                # Adding 1 to because of shifted geometric in scipy

                # Parameters
                n_states = 4
                label = {0: "ES", 1: "GD", 2: "NE", 3: "GA"}

                reads = O - 1
                reads_nz = sorted(reads[reads != 0])
                size = len(reads_nz)
                mean_r = numpy.average(reads_nz[: int(0.95 * size)])
                mu = numpy.array([1 / 0.99, 0.01 * mean_r + 2, mean_r, mean_r * 5.0])
                # mu = numpy.array([1/0.99, 0.1 * mean_r + 2,  mean_r, mean_r*5.0])
                L = 1.0 / mu
                B = []  # Emission Probability Distributions
                for i in range(n_states):
                    B.append(scipy.stats.geom(L[i]).pmf)

                pins = self.calculate_pins(O - 1)
                pins_obs = sum([1 for rd in O if rd >= 2]) / float(len(O))
                pnon = 1.0 - pins
                pnon_obs = 1.0 - pins_obs

                for r in range(100):
                    if pnon ** r < 0.01:
                        break

                A = numpy.zeros((n_states, n_states))
                a = math.log1p(-B[int(n_states / 2)](1) ** r)
                b = r * math.log(B[int(n_states / 2)](1)) + math.log(
                    1.0 / 3
                )  # change to n_states-1?
                for i in range(n_states):
                    A[i] = [b] * n_states
                    A[i][i] = a

                PI = numpy.zeros(n_states)  # Initial state distribution
                PI[0] = 0.7
                PI[1:] = 0.3 / (n_states - 1)

                

                ###############
                ### VITERBI ###
                (Q_opt, delta, Q) = self.viterbi(A, B, PI, O)
                ###############

                ##################
                ### ALPHA PASS ###
                (log_prob_obs, alpha, C) = self.forward_procedure(numpy.exp(A), B, PI, O)
                ##################

                #################
                ### BETA PASS ###
                beta = self.backward_procedure(numpy.exp(A), B, PI, O, C)
                logging.log("Finished backward_procedure")
                #################

                T = len(O)
                total = 0
                state2count = dict.fromkeys(range(n_states), 0)
                for t in range(T):
                    state = Q_opt[t]
                    state2count[state] += 1
                    total += 1
                
                rows = []
                states = [int(Q_opt[t]) for t in range(T)]
                last_orf = ""
                for t in range(T):
                    percentage = (t/T)*100
                    progress_update(f"Generating output: {percentage:3.1f}%", percentage)
                    s_lab = label.get(states[t], "Unknown State")
                    gamma_t = (alpha[:, t] * beta[:, t]) / numpy.sum(alpha[:, t] * beta[:, t])
                    genes_at_site = hash.get(position[t], [""])
                    genestr = ""
                    if not (len(genes_at_site) == 1 and not genes_at_site[0]):
                        genestr = ",".join(
                            ["%s_(%s)" % (g, rv2info.get(g, "-")[0]) for g in genes_at_site]
                        )
                    
                    rows.append((
                        int(position[t]),
                        int(O[t]) - 1,
                        *gamma_t,
                        s_lab,
                        genestr,
                    ))
                
                # 
                # confidences
                # 
            # 
            # Write outputs
            # 
            if True:
                base_path = self.inputs.output_path
                if base_path: # path=None means write to STDOUT
                    output_path = misc.inject_path_extension(base_path, extension="sites")
                    genes_path = misc.inject_path_extension(base_path, extension="genes")
                
                transit_tools.write_result(
                    path=output_path, # path=None means write to STDOUT
                    file_kind=SitesFile.identifier,
                    rows=rows,
                    column_names=SitesFile.column_names,
                    extra_info=dict(
                        calculation_time=f"{(time.time() - start_time):0.1f}seconds",
                        analysis_type=Method.identifier,
                        files=dict(
                            combined_wig=Method.inputs.data_sources[0],
                            annotation_path=Method.inputs.annotation_path,
                        ),
                        stats=dict(
                            mean=float(numpy.average(reads_nz)),
                            median=float(numpy.median(reads_nz)),
                            pins_observed=float(pins_obs),
                            pins_estimated=float(pins),
                            run_length=int(r),
                            state_means=                    "   ".join(["%s: %8.4f"   % (label[i], mu[i]  ) for i in range(n_states)]),
                            self_transition_probabilities=  "   ".join(["%s: %2.4e"   % (label[i], A[i][i]) for i in range(n_states)]),
                            state_emmision_parameters_theta="   ".join(["%s: %1.4f"   % (label[i], L[i]   ) for i in range(n_states)]),
                            state_distributions=            "   ".join(["%s: %2.2f%%" % (label[i], state2count[i] * 100.0 / total) for i in range(n_states) ]),
                        ),
                        parameters=dict(
                            replicates=self.inputs.replicates,
                            normalization=self.inputs.normalization,
                            loess_correction=self.inputs.loess_correction,
                            n_terminus=self.inputs.n_terminus,
                            c_terminus=self.inputs.c_terminus,
                            annotation_path=self.inputs.annotation_path,
                            output_path=output_path,
                        ),


                    ),
                )
                logging.log(f"Finished HMM - Sites: {output_path}")
                results_area.add(output_path)
                
                # 
                # Gene Files
                # 
                if True:
                    logging.log("Creating HMM Genes Level Output")
                    
                    temp_obs = numpy.zeros((1, len(O)))
                    temp_obs[0, :] = O - 1
                    self.post_process_genes(temp_obs, position, states, genes_path)
                    
                    logging.log("Adding File: %s" % (genes_path))
                    results_area.add(genes_path)
                    logging.log("Finished HMM Method")

    def forward_procedure(self, A, B, PI, O):
        logging.log("Starting HMM forward_procedure")
        T = len(O)
        N = len(B)
        alpha = numpy.zeros((N, T))
        C = numpy.zeros(T)

        alpha[:, 0] = PI * [B[i](O[0]) for i in range(N)]

        C[0] = 1.0 / numpy.sum(alpha[:, 0])
        alpha[:, 0] = C[0] * alpha[:, 0]

        for t in range(1, T):
            # B[i](O[:,t])  =>  numpy.prod(B[i](O[:,t]))
            # b_o = numpy.array([numpy.prod(B[i](O[:,t])) for i in range(N)])
            b_o = [B[i](O[t]) for i in range(N)]

            alpha[:, t] = numpy.dot(alpha[:, t - 1], A) * b_o

            C[t] = numpy.nan_to_num(1.0 / numpy.sum(alpha[:, t]))
            alpha[:, t] = numpy.nan_to_num(alpha[:, t] * C[t])

            if numpy.sum(alpha[:, t]) == 0:
                alpha[:, t] = 0.0000000000001

            percentage = (100.0 * self.count / self.max_iterations)
            text = "Running HMM forward_procedure... %1.1f%% (process 2/3)" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1
        
        percentage = 100
        text = "Running HMM forward_procedure... %1.1f%% (process 2/3)" % percentage
        progress_update(text, percentage)
        
        log_prob_obs = -(numpy.sum(numpy.log(C)))
        return (log_prob_obs, alpha, C)

    def backward_procedure(self, A, B, PI, O, C=numpy.array([])):
        logging.log("Starting HMM backward_procedure")
        N = len(B)
        T = len(O)
        beta = numpy.zeros((N, T))

        beta[:, T - 1] = 1.0
        if C.any():
            beta[:, T - 1] = beta[:, T - 1] * C[T - 1]

        for t in range(T - 2, -1, -1):
            # B[i](O[:,t])  =>  numpy.prod(B[i](O[:,t]))
            # b_o = numpy.array([numpy.prod(B[i](O[:,t])) for i in range(N)])
            b_o = [B[i](O[t]) for i in range(N)]

            beta[:, t] = numpy.nan_to_num(numpy.dot(A, (b_o * beta[:, t + 1])))

            if sum(beta[:, t]) == 0:
                beta[:, t] = 0.0000000000001

            if C.any():
                beta[:, t] = beta[:, t] * C[t]
            
            percentage = (100.0 * self.count / self.max_iterations)
            text = "Running HMM backward_procedure ... %1.1f%% (process 3/3)" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1

        return beta

    def viterbi(self, A, B, PI, O):
        logging.log("Starting HMM viterbi")
        N = len(B)
        T = len(O)
        delta = numpy.zeros((N, T))

        b_o = [B[i](O[0]) for i in range(N)]
        delta[:, 0] = numpy.log(PI) + numpy.log(b_o)

        Q = numpy.zeros((N, T), dtype=int)

        numpy.seterr(divide="ignore")
        for t in range(1, T):
            b_o = [B[i](O[t]) for i in range(N)]
            # nus = delta[:, t-1] + numpy.log(A)
            nus = delta[:, t - 1] + A
            delta[:, t] = nus.max(1) + numpy.log(b_o)
            Q[:, t] = nus.argmax(1)
            
            percentage = (100.0 * self.count / self.max_iterations)
            text = "Running HMM viterbi... %5.1f%% (process 1/3)" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1

        Q_opt = [int(numpy.argmax(delta[:, T - 1]))]
        for t in range(T - 2, -1, -1):
            Q_opt.insert(0, Q[Q_opt[0], t + 1])

            percentage = (100.0 * self.count / self.max_iterations)
            text = "Running HMM viterbi... %5.1f%% (process 1/3)" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1

        numpy.seterr(divide="warn")
        percentage = (100.0 * self.count / self.max_iterations)
        text = "Running HMM viterbi... %5.1f%% (process 1/3)" % percentage
        if self.count % 1000 == 0:
            progress_update(text, percentage)

        return (Q_opt, delta, Q)

    def calculate_pins(self, reads):
        logging.log("Calculating pins")
        non_ess_reads = []
        temp = []
        for rd in reads:

            if rd >= 1:
                if len(temp) < 10:
                    non_ess_reads.extend(temp)
                non_ess_reads.append(rd)
                temp = []
            else:
                temp.append(rd)

        return sum([1 for rd in non_ess_reads if rd >= 1]) / float(len(non_ess_reads))

    def post_process_genes(self, data, position, states, output_path):
        with gui_tools.nice_error_log:
            # 
            # organize data
            # 
            if True:
                pos2state = dict([(position[t], states[t]) for t in range(len(states))])
                theta = numpy.mean(data > 0)
                genes = tnseq_tools.Genes(
                    tnseq_tools.CombinedWig.PositionsAndReads((self.inputs.ctrl_read_counts, self.inputs.ctrl_positions)),
                    self.inputs.annotation_path,
                    data=data,
                    position=position,
                    ignore_codon=False,
                    n_terminus=self.inputs.n_terminus,
                    c_terminus=self.inputs.c_terminus,
                )

                num2label = {0: "ES", 1: "GD", 2: "NE", 3: "GA"}

                rows = []
                gene_calls = {
                    "ES":0,
                    "GD":0,
                    "NE":0,
                    "GA":0,
                    "N/A":0
                }
                for each_gene in genes:
                    reads_nz = [c for c in each_gene.reads.flatten() if c > 0]
                    avg_read_nz = 0
                    if len(reads_nz) > 0:
                        avg_read_nz = numpy.average(reads_nz)

                    # State
                    gene_states = [pos2state[p] for p in each_gene.position]
                    statedist = {}
                    for st in gene_states:
                        if st not in statedist:
                            statedist[st] = 0
                        statedist[st] += 1

                    # State gene_calls
                    n0 = statedist.get(0, 0)
                    n1 = statedist.get(1, 0)
                    n2 = statedist.get(2, 0)
                    n3 = statedist.get(3, 0)

                    if each_gene.n > 0:
                        # this was intended to call genes ES if have sufficiently long run, but n0 (#ES) not even consecutive
                        # E = tnseq_tools.expected_runs(each_gene.n,   1.0 - theta)
                        # V = tnseq_tools.variance_run(each_gene.n,   1.0 - theta)
                        if n0 == each_gene.n:
                            gene_category = "ES"
                        # elif n0 >= int(E+(3*math.sqrt(V))): gene_category = "ES"
                        else:
                            temp = max([(statedist.get(s, 0), s) for s in [0, 1, 2, 3]])[1]
                            gene_category = num2label[temp]
                    else:
                        E = 0.0
                        V = 0.0
                        gene_category = "N/A"
                    
                    
                    rows.append(
                        (
                            each_gene.orf,
                            each_gene.name,
                            each_gene.desc,
                            each_gene.n,
                            n0,
                            n1,
                            n2,
                            n3,
                            "%1.4f" % each_gene.theta(),
                            "%1.2f" % avg_read_nz,
                            gene_category,
                        )
                    )
                    if gene_category not in gene_calls:
                        gene_calls[gene_category] = 0
                    gene_calls[gene_category] += 1
                
                # 
                # add confidence values if needed
                # 
                column_names = GeneFile.column_names
                extra_data_info = {}
                if self.inputs.conf_on:
                    column_names = column_names + HmmConfidenceHelper.column_names
                    row_extensions, conf_comments = HmmConfidenceHelper.compute_row_extensions(rows)
                    from collections import Counter
                    assert len(rows)==len(row_extensions),"row_extensions should match number of rows"
                    for index, (each_row, each_extension) in enumerate(zip(rows, row_extensions)):
                        rows[index] = tuple(each_row) + tuple(each_extension)
                    conf_info = {}
                    for x in conf_comments: y = x.split(":"); conf_info[y[0]] = y[1]
                    extra_data_info["Confidence Summary"] = conf_info
                    
                    # TODO:
                        # add metrics: 
                        # - low-confidence count (flags count)
                        # - ambiguous count (flags count)
            # 
            # Write data
            # 
            transit_tools.write_result(
                path=output_path, # path=None means write to STDOUT
                file_kind=GeneFile.identifier,
                rows=rows,
                column_names=column_names,
                extra_info={
                    "Summary Of Gene Calls": gene_calls,
                    "Naming Reference": {
                        "ES": "essential",
                        "GD": "insertions cause growth-defect",
                        "NE": "non-essential",
                        "GA": "insertions confer growth-advantage",
                        "N/A": "not analyzed (genes with 0 TA sites)",
                    },
                    **extra_data_info,
                },
            )


import numpy
class HmmConfidenceHelper:
    states = ["ES","GD","NE","GA"]
    column_names = [
        "Mean",
        "Consistency",
        "ES Probability",
        "GD Probability",
        "NE Probability",
        "GA Probability",
        "Confidence",
        "Flag",
    ]

    def calc_probs(nzmean, mean_params, sat):
        probs = []
        for state in HmmConfidenceHelper.states:
            mean_means, std_means = mean_params[state]
            probs.append(
                0
                if mean_means < 0
                else scipy.stats.norm.pdf(sat * nzmean, loc=mean_means, scale=std_means)
            )
        
        def normalize(L):
            tot = float(sum(L))
            return [x / tot for x in L]
        
        return normalize(probs)
    
    def first_pass(rows):
        consistency_values = []
        mean_insertions_per_gene = {}
        nz_means_per_gene = {}
        calls_per_gene = {}
        mean_something_per_gene = {}
        ta_site_count_per_gene = {}
        for orf, gene_name, description, total_sites, es_count, gd_count, ne_count, ga_count, mean_insertions, mean_reads, state_call in rows:
            total_sites = int(total_sites)
            if total_sites == 0:
                continue
            
            mean_insertions_per_gene[orf] = float(mean_insertions)
            nz_means_per_gene[orf]        = float(mean_reads)
            mean_something_per_gene[orf]  = mean_insertions_per_gene[orf] * nz_means_per_gene[orf]
            ta_site_count_per_gene[orf]   = total_sites
            calls_per_gene[orf]           = state_call
            votes = [int(x) for x in [ es_count, gd_count, ne_count, ga_count ]]
            consistency = max(votes) / float(total_sites)
            consistency_values.append(consistency)

        return consistency_values, calls_per_gene, nz_means_per_gene, mean_insertions_per_gene, mean_something_per_gene, ta_site_count_per_gene
    
    def compute_row_extensions(
        rows,
        pseudo_count=0.01, # shrink range of saturation form [0,1] to [0.01,0.99] to prevent nan's from Beta
    ):
        consistency_values, calls_per_gene, nz_means_per_gene, mean_insertions_per_gene, mean_something_per_gene, ta_site_count_per_gene = HmmConfidenceHelper.first_pass(rows)
        
        output_comments = ["# HMM confidence info:"]
        output_comments.append(
            "# avg gene-level consistency of HMM states: %s" % (round(numpy.mean(consistency_values), 4))
        )

        orf_ids = [x[0] for x in rows]

        sat_params = {}
        nz_mean_params = {}
        mean_params = {}

        output_comments.append("# state posterior probability distributions:")
        for state in HmmConfidenceHelper.states:
            sub = [ orf for orf in orf_ids if calls_per_gene.get(orf) == state ]
            nzmeans = [nz_means_per_gene[orf] for orf in sub]
            sats = [mean_insertions_per_gene[orf] for orf in sub]
            means = [mean_something_per_gene[orf] for orf in sub]

            mean_sat = numpy.mean(sats)
            std_sat = numpy.std(sats)
            med_nz_means = numpy.median(nzmeans)
            iqr_nz_means = scipy.stats.iqr(nzmeans)
            mean_means = (
                -999 if len(means) == 0 else numpy.median(means)
            )  # -999 if there are no GA genes, for example
            std_means = max(
                1.0, 0.7314 * scipy.stats.iqr(means)
            )  # don't let stdev collapse to 0 for ES

            # model nz_mean with robust Normal distribution: use median and IQR
            sigma = 0.7413 * iqr_nz_means
            sigma = max(0.01, sigma)  # don't let it collapse to 0
            nz_mean_params[state] = (med_nz_means, sigma)

            # model saturation as Beta distribution, fit by method of moments:
            # https://real-statistics.com/distribution-fitting/method-of-moments/method-of-moments-beta-distribution/
            alpha = mean_sat * ((mean_sat * (1.0 - mean_sat) / (std_sat * std_sat) - 1.0))
            beta = alpha * (1.0 - mean_sat) / mean_sat
            sat_params[state] = (alpha, beta)
            mean_params[state] = (mean_means, std_means)

            output_comments.append(
                "#   Mean[%s]:   Norm(mean=%s,stdev=%s)"
                % (state, round(mean_params[state][0], 2), round(mean_params[state][1], 2))
            )

        # prob of each state is based Gaussian density for Mean count for each gene (combines Sat and nz_mean)

        # second pass...
        num_low_conf, num_ambig = 0, 0
        row_extensions = []
        for orf, gene_name, description, total_sites, es_count, gd_count, ne_count, ga_count, mean_insertions, mean_reads, state_call in rows:
            total_sites = int(total_sites)
            if total_sites == 0:
                row_extensions.append([""]*7) # blanks, to keep in register with rows (genes from HMM)
                continue

            votes = [int(x) for x in [es_count, gd_count, ne_count, ga_count]]
            consistency = max(votes) / float(total_sites)
            sat, nz_mean = float(mean_insertions), float(mean_reads)
            sat = max(pseudo_count, min(1.0 - pseudo_count, sat))
            probs = HmmConfidenceHelper.calc_probs(nz_mean, mean_params, sat)  # normalized
            conf = probs[HmmConfidenceHelper.states.index(state_call)]
            
            flag = ""

            if conf < 0.2:
                flag = "low-confidence"
            if conf >= 0.2 and conf != max(probs):
                flag = "ambiguous"

            if flag == "ambiguous":
                num_ambig += 1
            if flag == "low-confidence":
                num_low_conf += 1

            row_extension = (
                [round(sat * nz_mean, 1), round(consistency, 3)]
                + [round(x, 6) for x in probs]
            )
            row_extension += [round(conf, 4), flag]
            row_extensions.append(row_extension)

        output_comments.append("# num low-confidence genes: %s" % num_low_conf)
        output_comments.append("# num ambiguous genes: %s" % num_ambig)

        return row_extensions, output_comments # it would have been safer to output extensions appended onto rows

@transit_tools.ResultsFile
class SitesFile:
    identifier = "HMM_Sites"
    column_names = [
        "Location",
        "Read Count",
        "Probability ES",
        "Probability GD",
        "Probability NE",
        "Probability GA",
        "State",
        "Gene",
    ]
    
    @classmethod
    def can_load(cls, path):
        return transit_tools.file_starts_with(path, '#'+cls.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=self.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=self.identifier,
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[ "Probability ES" ]
                ).Show(),
            })
        )
        
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        #self.values_for_result_table.update({" ":""})

        parameters = self.extra_data.get("parameters",{})
        parameters_str = [str(key)+" : "+str(parameters[key]) for key in ["replicates","normalization", "loess_correction"]]
        self.values_for_result_table.update({"parameters": "; ".join(parameters_str) })

            
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

@transit_tools.ResultsFile
class GeneFile:
    identifier = "HMM_Genes"
    column_names = [
        "ORF",
        "Gene Name",
        "Description",
        "Total Sites",
        "ES Count",
        "GD Count",
        "NE Count",
        "GA Count",
        "Saturation",
        "NZmean",
        "State Call",
    ]
    
    @classmethod
    def can_load(cls, path):
        return transit_tools.file_starts_with(path, '#'+cls.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=self.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=self.identifier,
                    heading=self.comments_string,
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[
                        "ORF"
                    ],
                ).Show(),
            })
        )
        
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        
        # 
        # get summary stats
        #
        summary = self.extra_data.get("Summary Of Gene Calls", {})
        summary_str = [str(summary[key])+" "+str(key) for key in ["ES", "GD", "GA", "NE", "N/A"]] 
        self.values_for_result_table.update({"summary": "; ".join(summary_str) })
    
        

    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()


