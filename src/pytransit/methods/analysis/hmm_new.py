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
from pytransit.basics.lazy_dict import LazyDict

from pytransit.methods import analysis_base as base
from pytransit.tools.transit_tools import wx, pub, basename, HAS_R, FloatVector, DataFrame, StrVector, EOL
from pytransit.tools.tnseq_tools import Wig
from pytransit.tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools
from pytransit.universal_data import universal
from pytransit.components.parameter_panel import progress_update
from pytransit.components.spreadsheet import SpreadSheet
import pytransit.basics.csv as csv
import pytransit.components.file_display as file_display
import pytransit.components.samples_area as samples_area
import pytransit.components.results_area as results_area
import pytransit.basics.misc as misc
command_name = sys.argv[0]

# Checklist
    # DONE: basic data
    # DONE: create new inputs
    # DONE: translate the old define_panel
    # DONE: translate the old from_gui except for ctrl files (replace them with ctrl_read_counts and ctrl_positions)
    # DONE: change variable names to use new inputs
    # DONE: change from_args to set ctrl_read_counts, ctrl_positions
    # DONE: change Run() to use ctrl_read_counts, ctrl_positions
    # DONE: add all supporting methods
    # TODO: define file types 
    
@misc.singleton
class Analysis:
    short_name = "hmm"
    long_name = "HMM"
    short_desc = "Analysis of genomic regions using a Hidden Markov Model"
    long_desc = """Analysis of essentiality in the entire genome using a Hidden Markov Model. Capable of determining regions with different levels of essentiality representing Essential, Growth-Defect, Non-Essential and Growth-Advantage regions. Reference: DeJesus et al. (2013; BMC Bioinformatics)"""

    transposons = ["himar1"]
    inputs = LazyDict(
        data_sources=[],
        ctrl_read_counts=None,
        ctrl_positions=None,
        
        annotation_path=None,
        output_path=None,
        replicates="Mean",
        normalization=None,
        loess_correction=False,
        ignore_codon=True,
        n_terminus=0.0,
        c_terminus=0.0,
    )
    
    valid_cli_flags = [
        "-r",
        "-n",
        "-l",
        "-iN",
        "-iC",
    ]
    
    usage_string = f"""python3 {sys.argv[0]} hmm <comma-separated .wig files> <annotation .prot_table or GFF3> <output file>

        Optional Arguments:
            -r <string>     :=  How to handle replicates. Sum, Mean. Default: -r Mean
            -n <string>     :=  Normalization method. Default: -n TTR
            -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Off.
            -iN <float>     :=  Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 0
            -iC <float>     :=  Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 0
    """.replace("\n        ", "\n")
    
    
    wxobj = None
    panel = None
    
    def __init__(self, *args, **kwargs):
        self.instance = self.method = self.gui = self # for compatibility with older code/methods
        self.full_name        = f"[{self.short_name}]  -  {self.short_desc}"
        self.transposons_text = transit_tools.get_transposons_text(self.transposons)
    
    def __str__(self):
        return f"""
            Analysis Method:
                Short Name:  {self.short_name}
                Long Name:   {self.long_name}
                Short Desc:  {self.short_desc}
                Long Desc:   {self.long_desc}
                GUI:         {self.gui}
        """.replace('\n            ','\n').strip()
    
    def __repr__(self): return f"{self.inputs}"
    def __call__(self): return self

    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (self.panel, main_sizer):
            # 
            # parameter inputs
            # 
            self.value_getters = LazyDict()
            if True:
                self.value_getters.normalization          = panel_helpers.create_normalization_input(self.panel, main_sizer)
                self.value_getters.replicates             = panel_helpers.create_choice_input(self.panel, main_sizer, label="Replicates:", options=["Mean", "Sum"], tooltip_text="Determines how to handle replicates, and their read-counts. When using many replicates, using 'Mean' may be recommended over 'Sum'")
                self.value_getters.condition              = panel_helpers.create_condition_input(self.panel, main_sizer)
                self.value_getters.n_terminus             = panel_helpers.create_n_terminus_input(self.panel, main_sizer)
                self.value_getters.c_terminus             = panel_helpers.create_c_terminus_input(self.panel, main_sizer)
                self.value_getters.loess_correction       = panel_helpers.create_check_box_getter(self.panel, main_sizer, label_text="Correct for Genome Positional Bias", default_value=False, tooltip_text="Check to correct read-counts for possible regional biase using loess_correction. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
                
                panel_helpers.create_run_button(self.panel, main_sizer, from_gui_function=self.from_gui)
        
    @staticmethod
    def from_gui(frame):
        # 
        # get wig files
        # 
        combined_wig = universal.session_data.combined_wigs[0]
        Analysis.inputs.combined_wig = combined_wig.main_path
        Analysis.inputs.metadata     = combined_wig.metadata.path
        
        # 
        # get annotation
        # 
        Analysis.inputs.annotation_path = universal.session_data.annotation_path
        if not transit_tools.validate_annotation(Analysis.inputs.annotation_path):
            return None
        
        # TODO: validate transposons
        
        # 
        # setup custom inputs
        # 
        for each_key, each_getter in Analysis.instance.value_getters.items():
            try:
                Analysis.inputs[each_key] = each_getter()
            except Exception as error:
                raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
        
        # 
        # validate
        # 
        assert Analysis.inputs.condition != "[None]", "Please select a condition"
        
        # 
        # save result files
        # 
        Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name="hmm_output.dat",
            output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
        )
        if not Analysis.inputs.output_path:
            return None
        
        Analysis.inputs.data_sources = [ universal.session_data.combined_wigs[0].main_path ]
        
        # 
        # extract universal data
        # 
        Analysis.inputs.ctrl_read_counts, Analysis.inputs.ctrl_positions = transit_tools.gather_sample_data_for(conditions=[ Analysis.instance.value_getters.condition() ])
        
        return Analysis.instance

    @staticmethod
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Analysis.usage_string)
        console_tools.handle_unrecognized_flags(Analysis.valid_cli_flags, kwargs, Analysis.usage_string)

        ctrldata        = args[0].split(",")
        annotation_path = args[1]
        output_path     = args[2]
        
        replicates    = kwargs.get("r", Analysis.inputs.replicates)
        normalization = kwargs.get("n", Analysis.inputs.normalization)
        loess_correction         = kwargs.get("l", Analysis.inputs.loess_correction)
        ignore_codon  = True # TODO: not sure why this is hardcoded --Jeff
        n_terminus    = float(kwargs.get("iN", Analysis.inputs.n_terminus))
        c_terminus    = float(kwargs.get("iC", Analysis.inputs.c_terminus))
        
        (ctrl_read_counts, ctrl_positions) = transit_tools.get_validated_data(ctrldata)

        Analysis.inputs.update(dict(
            data_sources=ctrldata,
            ctrl_read_counts=ctrl_read_counts,
            ctrl_positions=ctrl_positions,
            annotation_path=annotation_path,
            output_path=output_path,
            replicates=replicates,
            normalization=normalization,
            loess_correction=loess_correction,
            ignore_codon=ignore_codon,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
        ))
        return Analysis.instance

    def Run(self):
        with gui_tools.nice_error_log:
            # 
            # Calculations
            # 
            if True:
                self.max_iterations = len(Analysis.inputs.ctrl_positions) * 4 + 1 # FIXME: I'm not sure this is right -- Jeff
                self.count = 1
                logging.log("Starting HMM Method")
                start_time = time.time()
                
                data, position = Analysis.inputs.ctrl_read_counts, Analysis.inputs.ctrl_positions
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
                    from pytransit.tools import stat_tools
                    for j in range(K):
                        data[j] = stat_tools.loess_correction(position, data[j])

                hash = transit_tools.get_pos_hash(self.inputs.annotation_path)
                rv2info = transit_tools.get_gene_info(self.inputs.annotation_path)

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
                    progress_update(f"Generating lines: {percentage:f3.2}%", percentage)
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
                        *gamma_t, # TODO: previously these numbers were formatted in a specific way
                        s_lab,
                        genestr,
                    ))
            # 
            # Write output
            # 
            if True:
                self.inputs.ctrl_read_counts = self.inputs.ctrl_read_counts.tolist()
                self.inputs.ctrl_positions   = self.inputs.ctrl_positions.tolist()
                extra_info = dict(
                        gui_or_cli=universal.interface,
                        cli_args=sys.argv,
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
                        parameters=self.inputs,
                    )
                transit_tools.write_result(
                    path=self.inputs.output_path, # path=None means write to STDOUT
                    file_kind=SitesFile.identifier,
                    rows=rows,
                    column_names=SitesFile.column_names,
                    extra_info=extra_info,
                )
                logging.log(f"Finished HMM - Sites: {self.inputs.output_path}")
                results_area.add(self.inputs.output_path)
            # 
            # Gene Files
            # 
            if True:
                logging.log("Creating HMM Genes Level Output")
                genes_path = (
                    ".".join(self.inputs.output_path.split(".")[:-1])
                    + "_genes."
                    + self.inputs.output_path.split(".")[-1]
                )

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
            text = "Running HMM forward_procedure... %1.1f%%" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1
            # print(t, O[:,t], alpha[:,t])
        
        percentage = 100
        text = "Running HMM forward_procedure... %1.1f%%" % percentage
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
            text = "Running HMM backward_procedure ... %1.1f%%" % percentage
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
            text = "Running HMM viterbi... %5.1f%%" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1

        Q_opt = [int(numpy.argmax(delta[:, T - 1]))]
        for t in range(T - 2, -1, -1):
            Q_opt.insert(0, Q[Q_opt[0], t + 1])

            percentage = (100.0 * self.count / self.max_iterations)
            text = "Running HMM viterbi... %5.1f%%" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1

        numpy.seterr(divide="warn")
        percentage = (100.0 * self.count / self.max_iterations)
        text = "Running HMM viterbi... %5.1f%%" % percentage
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
        logging.log("starting next step")
        with open(output_path, "w") as output:
            pos2state = dict([(position[t], states[t]) for t in range(len(states))])
            theta = numpy.mean(data > 0)
            G = tnseq_tools.Genes(
                tnseq_tools.CombinedWig.PositionsAndReads((self.inputs.ctrl_read_counts, self.inputs.ctrl_positions)),
                self.inputs.annotation_path,
                data=data,
                position=position,
                ignore_codon=False,
                n_terminus=self.inputs.n_terminus,
                c_terminus=self.inputs.c_terminus,
            )

            num2label = {0: "ES", 1: "GD", 2: "NE", 3: "GA"}
            output.write(f"#{GeneFile.identifier}\n")

            lines, counts = [], {}
            for gene in G:
                reads_nz = [c for c in gene.reads.flatten() if c > 0]
                avg_read_nz = 0
                if len(reads_nz) > 0:
                    avg_read_nz = numpy.average(reads_nz)

                # State
                genestates = [pos2state[p] for p in gene.position]
                statedist = {}
                for st in genestates:
                    if st not in statedist:
                        statedist[st] = 0
                    statedist[st] += 1

                # State counts
                n0 = statedist.get(0, 0)
                n1 = statedist.get(1, 0)
                n2 = statedist.get(2, 0)
                n3 = statedist.get(3, 0)

                if gene.n > 0:
                    # this was intended to call genes ES if have sufficiently long run, but n0 (#ES) not even consecutive
                    # E = tnseq_tools.expected_runs(gene.n,   1.0 - theta)
                    # V = tnseq_tools.variance_run(gene.n,   1.0 - theta)
                    if n0 == gene.n:
                        S = "ES"
                    # elif n0 >= int(E+(3*math.sqrt(V))): S = "ES"
                    else:
                        temp = max([(statedist.get(s, 0), s) for s in [0, 1, 2, 3]])[1]
                        S = num2label[temp]
                else:
                    E = 0.0
                    V = 0.0
                    S = "N/A"
                lines.append(
                    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.4f\t%1.2f\t%s\n"
                    % (
                        gene.orf,
                        gene.name,
                        gene.desc,
                        gene.n,
                        n0,
                        n1,
                        n2,
                        n3,
                        gene.theta(),
                        avg_read_nz,
                        S,
                    )
                )
                if S not in counts:
                    counts[S] = 0
                counts[S] += 1

            output.write("#command line: python3 %s\n" % (" ".join(sys.argv)))
            output.write("#summary of gene calls: ES=%s, GD=%s, NE=%s, GA=%s, N/A=%s\n"% tuple([counts.get(x, 0) for x in "ES GD NE GA N/A".split()]))
            output.write("#key: ES=essential, GD=insertions cause growth-defect, NE=non-essential, GA=insertions confer growth-advantage, N/A=not analyzed (genes with 0 TA sites)\n")
            output.write("#ORF\tgene\tannotation\tTAs\tES sites\tGD sites\tNE sites\tGA sites\tsaturation\tmean\tcall\n" )
            for line in lines:
                output.write(line)

@transit_tools.ResultsFile
class SitesFile:
    identifier = "HMM_Sites"
    column_names = [
        "Location",
        "Read Count",
        "Probability - ES",
        "Probability - GD",
        "Probability - NE",
        "Probability - GA",
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
            type=Analysis.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Analysis.identifier,
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[ "Adj. p-value", "p-value" ]
                ).Show(),
            })
        )
        
        self.column_names, self.rows, self.extra_data = tnseq_tools.read_results_file(self.path)
        self.values_for_result_table.update(self.extra_data.get("parameters", {}))
    
    def __str__(self):
        return f"""
            File for {Analysis.short_name}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

@transit_tools.ResultsFile
class GeneFile:
    identifier = "HMM_Genes"
    column_names = [
        "Orf",
        "Name",
        "Description",
        "Total Sites",
        "Num. ES",
        "Num. GD",
        "Num. NE",
        "Num. GA",
        "Avg. Insertions",
        "Avg. Reads",
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
            type=Analysis.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Analysis.identifier,
                    heading=self.comments,
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[
                        # HANDLE_THIS
                    ],
                ).Show(),
            })
        )
        
        # 
        # get column names
        # 
        comments, headers, rows = csv.read(self.path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
        self.comments = "\n".join(comments)
        
        # 
        # get rows
        #
        self.rows = []
        for each_row in rows:
            self.rows.append({
                each_column_name: each_cell
                    for each_column_name, each_cell in zip(self.column_names, each_row)
            })
        
        # 
        # get summary stats
        #
        self.values_for_result_table.update({
            # HANDLE_THIS (additional summary_info for results table)
            # examples:
                # f"Gene Count": len(self.rows),
                # f"Padj<{Analysis.significance_threshold}": len([
                #     1 for each in self.rows
                #         if each.get("Padj", 0) < Analysis.significance_threshold 
                # ]),
        })
    
    def __str__(self):
        return f"""
            File for {Analysis.short_name}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()


Method = GUI = Analysis