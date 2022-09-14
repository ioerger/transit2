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
from pytransit.components.parameter_panel import panel as parameter_panel
from pytransit.components.parameter_panel import panel, progress_update
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.panel_helpers import make_panel, create_run_button, create_normalization_input, create_reference_condition_input, create_include_condition_list_input, create_exclude_condition_list_input, create_n_terminus_input, create_c_terminus_input, create_pseudocount_input, create_winsorize_input, create_alpha_input, create_button, create_text_box_getter, create_button, create_check_box_getter, create_control_condition_input, create_experimental_condition_input, create_preview_loess_button, create_choice_input
import pytransit.basics.csv as csv
import pytransit.components.file_display as file_display
import pytransit.components.samples_area as samples_area
import pytransit.components.results_area as results_area
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
    

class Analysis:
    identifier  = "Hmm"
    short_name = "hmm - new"
    long_name = "HMM - new"
    short_desc = "Analysis of genomic regions using a Hidden Markov Model"
    long_desc = """Analysis of essentiality in the entire genome using a Hidden Markov Model. Capable of determining regions with different levels of essentiality representing Essential, Growth-Defect, Non-Essential and Growth-Advantage regions. Reference: DeJesus et al. (2013; BMC Bioinformatics)"""

    transposons = ["himar1"]
    columns_sites = [
        "Location",
        "Read Count",
        "Probability - ES",
        "Probability - GD",
        "Probability - NE",
        "Probability - GA",
        "State",
        "Gene",
    ]
    columns_genes = [
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
    
    inputs = LazyDict(
        data_sources=[],
        ctrl_read_counts=None,
        ctrl_positions=None,
        
        annotation_path=None,
        output_file=None,
        replicates="Mean",
        normalization=None,
        loess=False,
        ignore_codon=True,
        n_terminus=0.0,
        c_terminus=0.0,
    )
    
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

    def define_panel(self, _):
        self.panel = make_panel()
        
        # 
        # parameter inputs
        # 
        self.value_getters = LazyDict()
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        if True:
            self.value_getters.normalization          = create_normalization_input(self.panel, main_sizer)
            self.value_getters.replicates             = create_choice_input(self.panel, main_sizer, label="Replicates:", options=["Mean", "Sum"], tooltip_text="Determines how to handle replicates, and their read-counts. When using many replicates, using 'Mean' may be recommended over 'Sum'")
            self.value_getters.selected_wig_id        = create_choice_input(self.panel, main_sizer, label="Wig data:", options=[ each["id"] for each in samples_area.sample_table.rows ])
            self.value_getters.condition              = create_control_condition_input(self.panel, main_sizer)
            self.value_getters.n_terminus             = create_n_terminus_input(self.panel, main_sizer)
            self.value_getters.c_terminus             = create_c_terminus_input(self.panel, main_sizer)
            self.value_getters.loess                  = create_check_box_getter(self.panel, main_sizer, label_text="Correct for Genome Positional Bias", default_value=False, tooltip_text="Check to correct read-counts for possible regional biase using loess. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
            create_preview_loess_button(self.panel, main_sizer, wig_ids_getter=lambda *args,**kwargs: [ self.value_getters.selected_wig_id() ])
            
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
            # save result files
            # 
            Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
                default_file_name="hmm_output.dat",
                output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
            )
            if not Analysis.inputs.output_path:
                return None
            
            Analysis.inputs.output_file = open(Analysis.inputs.output_path, "w")
            Analysis.inputs.data_sources = [ universal.session_data.combined_wigs[0].main_path ]
            
            # 
            # extract universal data
            # 
            Analysis.inputs.ctrl_read_counts, Analysis.inputs.ctrl_positions = transit_tools.gather_sample_data_for(conditions=[ Analysis.instance.value_getters.condition() ])
            
            return Analysis.instance

    @classmethod
    def from_args(cls, rawargs):
        (args, kwargs) = transit_tools.clean_args(rawargs)

        ctrldata = args[0].split(",")
        Analysis.inputs.annotation_path = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        replicates = kwargs.get("r", "Mean")
        normalization = kwargs.get("n", "TTR")
        loess = kwargs.get("l", False)
        ignore_codon = True
        n_terminus = float(kwargs.get("iN", 0.0))
        c_terminus = float(kwargs.get("iC", 0.0))
        
        (ctrl_read_counts, ctrl_positions) = transit_tools.get_validated_data(ctrldata)

        cls.inputs.update(dict(
            data_sources=ctrldata,
            ctrl_read_counts=ctrl_read_counts,
            ctrl_positions=ctrl_positions,
            annotation_path=Analysis.inputs.annotation_path,
            output_file=output_file,
            replicates=replicates,
            normalization=normalization,
            loess=loess,
            ignore_codon=ignore_codon,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
        ))
        return Analysis.instance

    def Run(self):
        with gui_tools.nice_error_log:
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

            # Do loess
            if self.inputs.loess:
                logging.log("Performing loess Correction")
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
            Nstates = 4
            label = {0: "ES", 1: "GD", 2: "NE", 3: "GA"}

            reads = O - 1
            reads_nz = sorted(reads[reads != 0])
            size = len(reads_nz)
            mean_r = numpy.average(reads_nz[: int(0.95 * size)])
            mu = numpy.array([1 / 0.99, 0.01 * mean_r + 2, mean_r, mean_r * 5.0])
            # mu = numpy.array([1/0.99, 0.1 * mean_r + 2,  mean_r, mean_r*5.0])
            L = 1.0 / mu
            B = []  # Emission Probability Distributions
            for i in range(Nstates):
                B.append(scipy.stats.geom(L[i]).pmf)

            pins = self.calculate_pins(O - 1)
            pins_obs = sum([1 for rd in O if rd >= 2]) / float(len(O))
            pnon = 1.0 - pins
            pnon_obs = 1.0 - pins_obs

            for r in range(100):
                if pnon ** r < 0.01:
                    break

            A = numpy.zeros((Nstates, Nstates))
            a = math.log1p(-B[int(Nstates / 2)](1) ** r)
            b = r * math.log(B[int(Nstates / 2)](1)) + math.log(
                1.0 / 3
            )  # change to Nstates-1?
            for i in range(Nstates):
                A[i] = [b] * Nstates
                A[i][i] = a

            PI = numpy.zeros(Nstates)  # Initial state distribution
            PI[0] = 0.7
            PI[1:] = 0.3 / (Nstates - 1)

            

            ###############
            ### VITERBI ###
            (Q_opt, delta, Q) = self.viterbi(A, B, PI, O)
            ###############

            ##################
            ### ALPHA PASS ###
            (log_Prob_Obs, alpha, C) = self.forward_procedure(numpy.exp(A), B, PI, O)
            ##################

            #################
            ### BETA PASS ###
            beta = self.backward_procedure(numpy.exp(A), B, PI, O, C)
            #################

            T = len(O)
            total = 0
            state2count = dict.fromkeys(range(Nstates), 0)
            for t in range(T):
                state = Q_opt[t]
                state2count[state] += 1
                total += 1

            self.output = self.inputs.output_file
            self.output.write("#HMM - Sites\n")
            self.output.write("# Tn-HMM\n")

            if self.wxobj:
                members = sorted(
                    [
                        attr
                        for attr in dir(self)
                        if not callable(getattr(self, attr)) and not attr.startswith("__")
                    ]
                )
                memberstr = ""
                for m in members:
                    memberstr += "%s = %s, " % (m, getattr(self, m))
                self.output.write(
                    "#GUI with: data_sources=%s, annotation=%s, output=%s\n"
                    % (
                        ",".join(self.inputs.data_sources).encode("utf-8"),
                        self.inputs.annotation_path.encode("utf-8"),
                        self.output.name.encode("utf-8"),
                    )
                )
            else:
                self.output.write("#Console: python3 %s\n" % " ".join(sys.argv))

            self.output.write("# \n")
            self.output.write("# Mean:\t%2.2f\n" % (numpy.average(reads_nz)))
            self.output.write("# Median:\t%2.2f\n" % numpy.median(reads_nz))
            self.output.write("# Normalization:\t%s\n" % self.inputs.normalization)
            self.output.write("# loess Correction:\t%s\n" % str(self.inputs.loess))
            self.output.write("# pins (obs):\t%f\n" % pins_obs)
            self.output.write("# pins (est):\t%f\n" % pins)
            self.output.write("# Run length (r):\t%d\n" % r)
            self.output.write("# State means:\n")
            self.output.write(
                "#    %s\n"
                % "   ".join(["%s: %8.4f" % (label[i], mu[i]) for i in range(Nstates)])
            )
            self.output.write("# Self-Transition Prob:\n")
            self.output.write(
                "#    %s\n"
                % "   ".join(["%s: %2.4e" % (label[i], A[i][i]) for i in range(Nstates)])
            )
            self.output.write("# State Emission Parameters (theta):\n")
            self.output.write(
                "#    %s\n"
                % "   ".join(["%s: %1.4f" % (label[i], L[i]) for i in range(Nstates)])
            )
            self.output.write("# State Distributions:")
            self.output.write(
                "#    %s\n"
                % "   ".join(
                    [
                        "%s: %2.2f%%" % (label[i], state2count[i] * 100.0 / total)
                        for i in range(Nstates)
                    ]
                )
            )

            states = [int(Q_opt[t]) for t in range(T)]
            last_orf = ""
            for t in range(T):
                s_lab = label.get(states[t], "Unknown State")
                gamma_t = (alpha[:, t] * beta[:, t]) / numpy.sum(alpha[:, t] * beta[:, t])
                genes_at_site = hash.get(position[t], [""])
                genestr = ""
                if not (len(genes_at_site) == 1 and not genes_at_site[0]):
                    genestr = ",".join(
                        ["%s_(%s)" % (g, rv2info.get(g, "-")[0]) for g in genes_at_site]
                    )

                self.output.write(
                    "%s\t%s\t%s\t%s\t%s\n"
                    % (
                        int(position[t]),
                        int(O[t]) - 1,
                        "\t".join(["%-9.2e" % g for g in gamma_t]),
                        s_lab,
                        genestr,
                    )
                )

            self.output.close()

            logging.log("")  # Printing empty line to flush stdout
            logging.log("Finished HMM - Sites Method")
            logging.log("Adding File: %s" % (self.output.name))
            results_area.add(self.output.name)

            # Gene Files
            logging.log("Creating HMM Genes Level Output")
            genes_path = (
                ".".join(self.output.name.split(".")[:-1])
                + "_genes."
                + self.output.name.split(".")[-1]
            )

            temp_obs = numpy.zeros((1, len(O)))
            temp_obs[0, :] = O - 1
            self.post_process_genes(temp_obs, position, states, genes_path)

            logging.log("Adding File: %s" % (genes_path))
            results_area.add(self.output.name)
            logging.log("Finished HMM Method")

    def forward_procedure(self, A, B, PI, O):
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
            text = "Running HMM Method... %1.1f%%" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1
            # print(t, O[:,t], alpha[:,t])

        log_Prob_Obs = -(numpy.sum(numpy.log(C)))
        return (log_Prob_Obs, alpha, C)

    def backward_procedure(self, A, B, PI, O, C=numpy.array([])):

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
            text = "Running HMM Method... %1.1f%%" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1

        return beta

    def viterbi(self, A, B, PI, O):
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
            text = "Running HMM Method... %5.1f%%" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1

        Q_opt = [int(numpy.argmax(delta[:, T - 1]))]
        for t in range(T - 2, -1, -1):
            Q_opt.insert(0, Q[Q_opt[0], t + 1])

            percentage = (100.0 * self.count / self.max_iterations)
            text = "Running HMM Method... %5.1f%%" % percentage
            if self.count % 1000 == 0:
                progress_update(text, percentage)
            self.count += 1

        numpy.seterr(divide="warn")
        percentage = (100.0 * self.count / self.max_iterations)
        text = "Running HMM Method... %5.1f%%" % percentage
        if self.count % 1000 == 0:
            progress_update(text, percentage)

        return (Q_opt, delta, Q)

    def calculate_pins(self, reads):
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
            output.write("#HMM - Genes\n")

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
            output.write(
                "#summary of gene calls: ES=%s, GD=%s, NE=%s, GA=%s, N/A=%s\n"
                % tuple([counts.get(x, 0) for x in "ES GD NE GA N/A".split()])
            )
            output.write(
                "#key: ES=essential, GD=insertions cause growth-defect, NE=non-essential, GA=insertions confer growth-advantage, N/A=not analyzed (genes with 0 TA sites)\n"
            )
            output.write(
                "#ORF\tgene\tannotation\tTAs\tES sites\tGD sites\tNE sites\tGA sites\tsaturation\tmean\tcall\n"
            )
            for line in lines:
                output.write(line)

@transit_tools.ResultsFile
class File(Analysis):
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Analysis.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Analysis.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title=Analysis.identifier,heading="",column_names=self.column_names,rows=self.rows, sort_by=[ "Adj. p-value", "p-value" ]).Show(),
                # "Display Heatmap": lambda *args: self.create_heatmap(infile=self.path, output_path=self.path+".heatmap.png"),
            })
        )
        
        # 
        # get column names
        # 
        comments, headers, rows = csv.read(self.path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
        if len(comments) == 0:
            raise Exception(f'''No comments in file, and I expected the last comment to be the column names, while to load Anova file "{self.path}"''')
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
    
    def create_heatmap(self, infile, output_path, topk=-1, qval=0.05, low_mean_filter=5):
        if not HAS_R:
            raise Exception(f'''Error: R and rpy2 (~= 3.0) required to run Heatmap''')
        headers = None
        data, hits = [], []
        number_of_conditions = -1

        with open(infile) as file:
            for line in file:
                w = line.rstrip().split("\t")
                if line[0] == "#" or (
                    "pval" in line and "padj" in line
                ):  # check for 'pval' for backwards compatibility
                    headers = w
                    continue  # keep last comment line as headers
                # assume first non-comment line is header
                if number_of_conditions == -1:
                    # ANOVA header line has names of conditions, organized as 3+2*number_of_conditions+3 (2 groups (means, LFCs) X number_of_conditions conditions)
                    number_of_conditions = int((len(w) - 6) / 2)
                    headers = headers[3 : 3 + number_of_conditions]
                    headers = [x.replace("Mean_", "") for x in headers]
                else:
                    means = [
                        float(x) for x in w[3 : 3 + number_of_conditions]
                    ]  # take just the columns of means
                    lfcs = [
                        float(x) for x in w[3 + number_of_conditions : 3 + number_of_conditions + number_of_conditions]
                    ]  # take just the columns of LFCs
                    each_qval = float(w[-2])
                    data.append((w, means, lfcs, each_qval))
        
        data.sort(key=lambda x: x[-1])
        hits, LFCs = [], []
        for k, (w, means, lfcs, each_qval) in enumerate(data):
            if (topk == -1 and each_qval < qval) or (
                topk != -1 and k < topk
            ):
                mm = round(numpy.mean(means), 1)
                if mm < low_mean_filter:
                    print("excluding %s/%s, mean(means)=%s" % (w[0], w[1], mm))
                else:
                    hits.append(w)
                    LFCs.append(lfcs)

        print("heatmap based on %s genes" % len(hits))
        genenames = ["%s/%s" % (w[0], w[1]) for w in hits]
        hash = {}
        headers = [h.replace("Mean_", "") for h in headers]
        for i, col in enumerate(headers):
            hash[col] = FloatVector([x[i] for x in LFCs])
        df = DataFrame(hash)
        transit_tools.r_heatmap_func(df, StrVector(genenames), output_path)
        
        # add it as a result
        results_area.add(output_path)
        gui_tools.show_image(output_path)

    def graph_volcano_plot(dataset_name, dataset_type, dataset_path):
        try:
            X = []
            Y = []
            header = []
            qval_list = []
            bad = []
            col_logFC = -6
            col_pval = -2
            col_qval = -1
            ii = 0
            with open(dataset_path) as file:
                for line in file:
                    if line.startswith("#"):
                        tmp = line.split("\t")
                        temp_col_logfc = [
                            i
                            for (i, x) in enumerate(tmp)
                            if "logfc" in x.lower()
                            or "log-fc" in x.lower()
                            or "log2fc" in x.lower()
                        ]
                        temp_col_pval = [
                            i
                            for (i, x) in enumerate(tmp)
                            if ("pval" in x.lower() or "p-val" in x.lower())
                            and "adj" not in x.lower()
                        ]
                        if temp_col_logfc:
                            col_logFC = temp_col_logfc[-1]
                        if temp_col_pval:
                            col_pval = temp_col_pval[-1]
                        continue

                    tmp = line.strip().split("\t")
                    try:
                        log10qval = -math.log(float(tmp[col_pval].strip()), 10)
                    except ValueError as e:
                        bad.append(ii)
                        log10qval = 0

                    log2FC = float(tmp[col_logFC])

                    qval_list.append( (float(tmp[col_qval]), float(tmp[col_pval].strip())) )
                    X.append(log2FC)
                    Y.append(log10qval)
                    ii += 1
            count = 0
            threshold = 0.00001
            backup_thresh = 0.00001
            qval_list.sort()
            for (q, p) in qval_list:
                backup_thresh = p
                if q > 0.05:
                    break
                threshold = p
                count += 1

            if threshold == 0:
                threshold = backup_thresh
            for ii in bad:
                Y[ii] = max(Y)
            plt.plot(X, Y, "bo")
            plt.axhline(
                -math.log(threshold, 10), color="r", linestyle="dashed", linewidth=3
            )
            plt.xlabel("Log Fold Change (base 2)")
            plt.ylabel("-Log p-value (base 10)")
            plt.suptitle("Resampling - Volcano plot")
            plt.title("Adjusted threshold (red line): %1.8f" % threshold)
            plt.show()

        except Exception as e:
            print("Error occurred creating plot:", str(e))

Method = GUI = Analysis
Analysis() # make sure there's one instance