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
import statsmodels.stats.multitest
from super_map import LazyDict

from pytransit.analysis import base
from pytransit.transit_tools import wx, pub, basename, HAS_R, FloatVector, DataFrame, StrVector, EOL
import pytransit
import pytransit.gui_tools as gui_tools
import pytransit.file_display as file_display
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools
import pytransit.basics.csv as csv
import pytransit.components.results_area as results_area
from pytransit.core_data import universal
from pytransit.components.parameter_panel import panel as parameter_panel
from pytransit.components.parameter_panel import panel, progress_update
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.panel_helpers import make_panel, create_run_button, create_normalization_input, create_reference_condition_input, create_include_condition_list_input, create_exclude_condition_list_input, create_n_terminus_input, create_c_terminus_input, create_pseudocount_input, create_winsorize_input, create_alpha_input, create_button, create_text_box_getter, create_button, create_check_box_getter
command_name = sys.argv[0]

class Analysis:
    identifier  = "#Resampling"
    short_name = "resampling - new"
    long_name = "Resampling (Permutation test)"
    short_desc = "Resampling test of conditional essentiality between two conditions"
    long_desc = """Method for determining conditional essentiality based on resampling (i.e. permutation test). Identifies significant changes in mean read-counts for each gene after normalization."""

    transposons = ["himar1", "tn5"]
    columns = [
        "Orf",
        "Name",
        "Desc",
        "Sites",
        "Mean Ctrl",
        "Mean Exp",
        "log2FC",
        "Sum Ctrl",
        "Sum Exp",
        "Delta Mean",
        "p-value",
        "Adj. p-value",
    ]
    
    inputs = LazyDict(
        # Fix these
        ctrldata=None,
        expdata=None,
        
        # 
        annotation_path=None,
        output_file=None,
        normalization="TTR",
        samples=10000,
        adaptive=False,
        do_histogram=False,
        include_zeros=False,
        pseudocount=1,
        replicates="Sum",
        LOESS=False,
        ignore_codon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        ctrl_lib_str="",
        exp_lib_str="",
        winz=False,
        Z=False,
        diff_strains=False,
        annotation_path_exp="",
        combined_wig_params=None,
    )
    
    valid_cli_flags = [
        "-n",
        "--include-conditions",
        "--exclude-conditions",
        "--ref",
        "-iN",
        "-iC",
        "-PC",
        "-alpha",
        "-winz",
    ]
    usage_string = f"""
        python3 {sys.argv[0]} resampling <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        ---
        OR
        ---
        python3 {sys.argv[0]} resampling -c <combined wig file> <samples_metadata file> <ctrl condition name> <exp condition name> <annotation .prot_table> <output file> [Optional Arguments]
        NB: The ctrl and exp condition names should match Condition names in samples_metadata file.

        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -n <string>     :=  Normalization method. Default: -n TTR
        -h              :=  Output histogram of the permutations for each gene. Default: Turned Off.
        -a              :=  Perform adaptive resampling. Default: Turned Off.
        -ez             :=  Exclude rows with zero across conditions. Default: Turned off
                            (i.e. include rows with zeros).
        -PC <float>     :=  Pseudocounts used in calculating LFC. (default: 1)
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias.
                            Default: Turned Off.
        -iN <int>       :=  Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 0
        -iC <int>       :=  Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 0
        --ctrl_lib      :=  String of letters representing library of control files in order
                            e.g. 'AABB'. Default empty. Letters used must also be used in --exp_lib
                            If non-empty, resampling will limit permutations to within-libraries.

        --exp_lib       :=  String of letters representing library of experimental files in order
                            e.g. 'ABAB'. Default empty. Letters used must also be used in --ctrl_lib
                            If non-empty, resampling will limit permutations to within-libraries.
        -winz           :=  winsorize insertion counts for each gene in each condition 
                            (replace max cnt in each gene with 2nd highest; helps mitigate effect of outliers)
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
            sample_getter          = create_text_box_getter(self.panel, main_sizer, label_text="Samples", default_value="10000", tooltip_text="Number of samples to take when estimating the resampling histogram. More samples give more accurate estimates of the p-values at the cost of computation time.")
            
            self.value_getters.samples                = lambda *args: int(sample_getter(*args))
            self.value_getters.n_terminus             = create_n_terminus_input(self.panel, main_sizer)
            self.value_getters.c_terminus             = create_c_terminus_input(self.panel, main_sizer)
            self.value_getters.pseudocount            = create_pseudocount_input(self.panel, main_sizer)
            self.value_getters.normalization          = create_normalization_input(self.panel, main_sizer)
            self.value_getters.genome_positional_bias = create_check_box_getter(self.panel, main_sizer, label_text="Correct for Genome Positional Bias", default_value=False, tooltip_text="Check to correct read-counts for possible regional biase using LOESS. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
            # FIXME: LOESS button create_button(panel, sizer, *, label)
            self.value_getters.adaptive               = create_check_box_getter(self.panel, main_sizer, label_text="Adaptive Resampling (Faster)", default_value=False, tooltip_text="Dynamically stops permutations early if it is unlikely the ORF will be significant given the results so far. Improves performance, though p-value calculations for genes that are not differentially essential will be less accurate.")
            self.value_getters.do_histogram            = create_check_box_getter(self.panel, main_sizer, label_text="Generate Resampling Histograms", default_value=False, tooltip_text="Creates .png images with the resampling histogram for each of the ORFs. Histogram images are created in a folder with the same name as the output file.")
            self.value_getters.include_zeros           = create_check_box_getter(self.panel, main_sizer, label_text="Include sites with all zeros", default_value=False, tooltip_text="Includes sites that are empty (zero) across all datasets. Unchecking this may be useful for tn5 datasets, where all nucleotides are possible insertion sites and will have a large number of empty sites (significantly slowing down computation and affecting estimates).")
            
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
            wig_group = universal.session_data.wig_groups[0]
            Analysis.inputs.combined_wig = wig_group.cwig.path
            Analysis.inputs.metadata     = wig_group.metadata.path
            
            # 
            # get annotation
            # 
            Analysis.inputs.annotation_path = universal.session_data.annotation_path
            # FIXME: enable this once I get a valid annotation file example
            # if not transit_tools.validate_annotation(Analysis.inputs.annotation):
            #     return None
            
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
                default_file_name="output.dat",
                output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
            )
            if not Analysis.inputs.output_path:
                return None

            return Analysis.instance

    @classmethod
    def from_args(self, rawargs):
        (args, kwargs) = transit_tools.clean_args(rawargs)

        isCombinedWig = True if kwargs.get("c", False) else False
        combined_wig_params = None
        if isCombinedWig:
            if len(args) != 5:
                print("Error: Incorrect number of args. See usage")
                print(self.usage_string)
                sys.exit(0)
            combined_wig_params = {
                "combined_wig": kwargs.get("c"),
                "samples_metadata": args[0],
                "conditions": [args[1].lower(), args[2].lower()],
            }
            annot_paths = args[3].split(",")
            # to show contrasted conditions for combined_wigs in output header
            ctrldata = [combined_wig_params["conditions"][0]]
            expdata = [combined_wig_params["conditions"][1]]
            output_path = args[4]
        else:
            if len(args) != 4:
                print("Error: Incorrect number of args. See usage")
                print(self.usage_string)
                sys.exit(0)
            ctrldata = args[0].split(",")
            expdata = args[1].split(",")
            annot_paths = args[2].split(",")
            output_path = args[3]
        annotation_path = annot_paths[0]
        diff_strains = False
        annotation_pathExp = ""
        if len(annot_paths) == 2:
            annotation_pathExp = annot_paths[1]
            diff_strains = True
        if diff_strains and isCombinedWig:
            print("Error: Cannot have combined wig and different annotation files.")
            sys.exit(0)
        winz = True if "winz" in kwargs else False

        output_file = open(output_path, "w")

        # check for unrecognized flags
        flags = (
            "-c -s -n -h -a -ez -PC -l -iN -iC --ctrl_lib --exp_lib -Z -winz".split()
        )
        for arg in rawargs:
            if arg[0] == "-" and arg not in flags:
                self.transit_error("flag unrecognized: %s" % arg)
                print(self.usage_string)
                sys.exit(0)

        normalization = kwargs.get("n", "TTR")
        samples = int(kwargs.get("s", 10000))
        adaptive = kwargs.get("a", False)
        do_histogram = kwargs.get("h", False)
        replicates = kwargs.get("r", "Sum")
        excludeZeros = kwargs.get("ez", False)
        include_zeros = not excludeZeros
        pseudocount = float(
            kwargs.get("PC", 1.0)
        )  # use -PC (new semantics: for LFCs) instead of -pc (old semantics: fake counts)

        Z = True if "Z" in kwargs else False

        LOESS = kwargs.get("l", False)
        ignore_codon = True

        n_terminus = float(kwargs.get("iN", 0.00))  # integer interpreted as percentage
        c_terminus = float(kwargs.get("iC", 0.00))
        ctrl_lib_str = kwargs.get("-ctrl_lib", "")
        exp_lib_str = kwargs.get("-exp_lib", "")
        
        self.inputs.update(dict(
            ctrldata=ctrldata,
            expdata=expdata,
            output_file=output_file,
            normalization=normalization,
            samples=samples,
            adaptive=adaptive,
            do_histogram=do_histogram,
            include_zeros=include_zeros,
            pseudocount=pseudocount,
            replicates=replicates,
            LOESS=LOESS,
            ignore_codon=ignore_codon,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            ctrl_lib_str=ctrl_lib_str,
            exp_lib_str=exp_lib_str,
            winz=winz,
            Z=Z,
            diff_strains=diff_strains,
            annotation_path_exp=annotation_path_exp,
            combined_wig_params=combined_wig_params,
        ))
        return self

    def Run(self):
        if self.do_histogram:
            try: import matplotlib.pyplot as plt
            except:
                print("Error: cannot do histograms")
                self.do_histogram = False

        transit_tools.log("Starting resampling Method")
        start_time = time.time()
        if self.winz:
            transit_tools.log("Winsorizing insertion counts")

        histPath = ""
        if self.do_histogram:
            histPath = os.path.join(
                os.path.dirname(self.output.name),
                transit_tools.fetch_name(self.output.name) + "_histograms",
            )
            if not os.path.isdir(histPath):
                os.makedirs(histPath)

        # Get orf data
        transit_tools.log("Getting Data")
        if self.diff_strains:
            transit_tools.log("Multiple annotation files found")
            transit_tools.log(
                "Mapping ctrl data to {0}, exp data to {1}".format(
                    self.annotation_path, self.annotation_path_exp
                )
            )

        if self.combined_wig_params:
            (position, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(
                self.combined_wig_params["combined_wig"]
            )
            conditionsByFile, _, _, _ = tnseq_tools.read_samples_metadata(
                self.combined_wig_params["samples_metadata"]
            )
            conditions = self.wigs_to_conditions(conditionsByFile, filenamesInCombWig)
            data, conditions = self.filter_wigs_by_conditions(
                data, conditions, self.combined_wig_params["conditions"]
            )
            data_ctrl = numpy.array(
                [
                    d
                    for i, d in enumerate(data)
                    if conditions[i].lower() == self.combined_wig_params["conditions"][0]
                ]
            )
            data_exp = numpy.array(
                [
                    d
                    for i, d in enumerate(data)
                    if conditions[i].lower() == self.combined_wig_params["conditions"][1]
                ]
            )
            position_ctrl, position_exp = position, position
        else:
            (data_ctrl, position_ctrl) = transit_tools.get_validated_data(
                self.ctrldata, wxobj=self.wxobj
            )
            (data_exp, position_exp) = transit_tools.get_validated_data(
                self.expdata, wxobj=self.wxobj
            )
        (K_ctrl, N_ctrl) = data_ctrl.shape
        (K_exp, N_exp) = data_exp.shape

        if not self.diff_strains and (N_ctrl != N_exp):
            self.transit_error(
                "Error: Ctrl and Exp wig files don't have the same number of sites."
            )
            self.transit_error("Make sure all .wig files come from the same strain.")
            return
        # (data, position) = transit_tools.get_validated_data(self.ctrldata+self.expdata, wxobj=self.wxobj)

        transit_tools.log("Preprocessing Ctrl data...")
        data_ctrl = self.preprocess_data(position_ctrl, data_ctrl)

        transit_tools.log("Preprocessing Exp data...")
        data_exp = self.preprocess_data(position_exp, data_exp)

        G_ctrl = tnseq_tools.Genes(
            self.ctrldata,
            self.annotation_path,
            ignore_codon=self.ignore_codon,
            n_terminus=self.n_terminus,
            c_terminus=self.c_terminus,
            data=data_ctrl,
            position=position_ctrl,
        )
        G_exp = tnseq_tools.Genes(
            self.expdata,
            self.annotation_path_exp,
            ignore_codon=self.ignore_codon,
            n_terminus=self.n_terminus,
            c_terminus=self.c_terminus,
            data=data_exp,
            position=position_exp,
        )

        doLibraryResampling = False
        # If library string not empty
        if self.ctrl_lib_str or self.exp_lib_str:
            letters_ctrl = set(self.ctrl_lib_str)
            letters_exp = set(self.exp_lib_str)

            # Check if using exactly 1 letters; i.e. no different libraries
            if len(letters_ctrl) == 1 and letters_exp == 1:
                pass
            # If using more than one letter, then check no differences in set
            else:
                lib_diff = letters_ctrl ^ letters_exp
                # Check that their differences
                if not lib_diff:
                    doLibraryResampling = True
                else:
                    transit_tools.transit_error(
                        "Error: Library Strings (Ctrl = %s, Exp = %s) do not use the same letters. Make sure every letter / library is represented in both Control and Experimental Conditions. Proceeding with resampling assuming all datasets belong to the same library."
                        % (self.ctrl_lib_str, self.exp_lib_str)
                    )
                    self.ctrl_lib_str = ""
                    self.exp_lib_str = ""

        (data, qval) = self.run_resampling(G_ctrl, G_exp, doLibraryResampling, histPath)
        self.write_output(data, qval, start_time)

        self.finish()
        transit_tools.log("Finished resampling Method")

    def preprocess_data(self, position, data):
        (K, N) = data.shape

        if self.normalization != "nonorm":
            transit_tools.log("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(
                data,
                self.normalization,
                self.ctrldata + self.expdata,
                self.annotation_path,
            )

        if self.LOESS:
            transit_tools.log("Performing LOESS Correction")
            for j in range(K):
                data[j] = stat_tools.loess_correction(position, data[j])

        return data

    def wigs_to_conditions(self, conditionsByFile, filenamesInCombWig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditionsByFile.get(f, None) for f in filenamesInCombWig]

    def filter_wigs_by_conditions(self, data, conditions, included_conditions):
        """
            Filters conditions from wig to ctrl, exp conditions only
            ([[Wigdata]], [ConditionCtrl, ConditionExp]) -> Tuple([[Wigdata]], [Condition])
        """
        d_filtered, cond_filtered = [], []
        if len(included_conditions) != 2:
            self.transit_error("Only 2 conditions expected", included_conditions)
            sys.exit(0)
        for i, c in enumerate(conditions):
            if c.lower() in included_conditions:
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])

        return (numpy.array(d_filtered), numpy.array(cond_filtered))

    def write_output(self, data, qval, start_time):

        self.output.write("#Resampling\n")
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
                "#GUI with: norm=%s, samples=%s, pseudocounts=%1.2f, adaptive=%s, histogram=%s, include_zeros=%s, output=%s\n"
                % (
                    self.normalization,
                    self.samples,
                    self.pseudocount,
                    self.adaptive,
                    self.do_histogram,
                    self.include_zeros,
                    self.output.name.encode("utf-8"),
                )
            )
        else:
            self.output.write("#Console: python3 %s\n" % " ".join(sys.argv))
        self.output.write(
            "#Parameters: samples=%s, norm=%s, histograms=%s, adaptive=%s, excludeZeros=%s, pseudocounts=%s, LOESS=%s, trim_Nterm=%s, trim_Cterm=%s\n"
            % (
                self.samples,
                self.normalization,
                self.do_histogram,
                self.adaptive,
                not self.include_zeros,
                self.pseudocount,
                self.LOESS,
                self.n_terminus,
                self.c_terminus,
            )
        )
        self.output.write(
            "#Control Data: %s\n" % (",".join(self.ctrldata).encode("utf-8"))
        )
        self.output.write(
            "#Experimental Data: %s\n" % (",".join(self.expdata).encode("utf-8"))
        )
        self.output.write(
            "#Annotation path: %s %s\n"
            % (
                self.annotation_path.encode("utf-8"),
                self.annotation_path_exp.encode("utf-8") if self.diff_strains else "",
            )
        )
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        # Z = True # include Z-score column in resampling output?
        global columns  # consider redefining columns above (for GUI)
        if self.Z == True:
            columns = [
                "Orf",
                "Name",
                "Desc",
                "Sites",
                "Mean Ctrl",
                "Mean Exp",
                "log2FC",
                "Sum Ctrl",
                "Sum Exp",
                "Delta Mean",
                "p-value",
                "Z-score",
                "Adj. p-value",
            ]
        self.output.write("#%s\n" % "\t".join(columns))

        for i, row in enumerate(data):
            (
                orf,
                name,
                desc,
                n,
                mean1,
                mean2,
                sum1,
                sum2,
                test_obs,
                log2FC,
                pval_2tail,
            ) = row
            if self.Z == True:
                p = pval_2tail / 2  # convert from 2-sided back to 1-sided
                if p == 0:
                    p = 1e-5  # or 1 level deeper the num of iterations of resampling, which is 1e-4=1/10000, by default
                if p == 1:
                    p = 1 - 1e-5
                z = scipy.stats.norm.ppf(p)
                if log2FC > 0:
                    z *= -1
                self.output.write(
                    "%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.1f\t%1.2f\t%1.1f\t%1.5f\t%0.2f\t%1.5f\n"
                    % (
                        orf,
                        name,
                        desc,
                        n,
                        mean1,
                        mean2,
                        log2FC,
                        sum1,
                        sum2,
                        test_obs,
                        pval_2tail,
                        z,
                        qval[i],
                    )
                )
            else:
                self.output.write(
                    "%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.1f\t%1.2f\t%1.1f\t%1.5f\t%1.5f\n"
                    % (
                        orf,
                        name,
                        desc,
                        n,
                        mean1,
                        mean2,
                        log2FC,
                        sum1,
                        sum2,
                        test_obs,
                        pval_2tail,
                        qval[i],
                    )
                )
        self.output.close()

        transit_tools.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)

    def winsorize_resampling(self, counts):
        # input is insertion counts for gene as pre-flattened numpy array
        counts = counts.tolist()
        if len(counts) < 3:
            return counts
        s = sorted(counts, reverse=True)
        if s[1] == 0:
            return counts  # don't do anything if there is only 1 non-zero value
        c2 = [s[1] if x == s[0] else x for x in counts]
        return numpy.array(c2)

    def run_resampling(
        self, G_ctrl, G_exp=None, doLibraryResampling=False, histPath=""
    ):
        data = []
        N = len(G_ctrl)
        count = 0
        self.progress_range(N)

        for gene in G_ctrl:
            if gene.orf not in G_exp:
                if self.diff_strains:
                    continue
                else:
                    self.transit_error(
                        "Error: Gene in ctrl data not present in exp data"
                    )
                    self.transit_error(
                        "Make sure all .wig files come from the same strain."
                    )
                    return ([], [])

            gene_exp = G_exp[gene.orf]
            count += 1

            if not self.diff_strains and gene.n != gene_exp.n:
                self.transit_error(
                    "Error: No. of TA sites in Exp and Ctrl data are different"
                )
                self.transit_error(
                    "Make sure all .wig files come from the same strain."
                )
                return ([], [])

            if (gene.k == 0 and gene_exp.k == 0) or gene.n == 0 or gene_exp.n == 0:
                (
                    test_obs,
                    mean1,
                    mean2,
                    log2FC,
                    pval_ltail,
                    pval_utail,
                    pval_2tail,
                    testlist,
                    data1,
                    data2,
                ) = (0, 0, 0, 0, 1.00, 1.00, 1.00, [], [0], [0])
            else:
                if not self.include_zeros:
                    ii_ctrl = numpy.sum(gene.reads, 0) > 0
                    ii_exp = numpy.sum(gene_exp.reads, 0) > 0
                else:
                    ii_ctrl = numpy.ones(gene.n) == 1
                    ii_exp = numpy.ones(gene_exp.n) == 1

                # data1 = gene.reads[:,ii_ctrl].flatten() + self.pseudocount # we used to have an option to add pseudocounts to each observation, like this
                data1 = gene.reads[:, ii_ctrl].flatten()
                data2 = gene_exp.reads[:, ii_exp].flatten()
                if self.winz:
                    data1 = self.winsorize_resampling(data1)
                    data2 = self.winsorize_resampling(data2)

                if doLibraryResampling:
                    (
                        test_obs,
                        mean1,
                        mean2,
                        log2FC,
                        pval_ltail,
                        pval_utail,
                        pval_2tail,
                        testlist,
                    ) = stat_tools.resampling(
                        data1,
                        data2,
                        S=self.samples,
                        testFunc=stat_tools.F_mean_diff_dict,
                        permFunc=stat_tools.F_shuffle_dict_libraries,
                        adaptive=self.adaptive,
                        lib_str1=self.ctrl_lib_str,
                        lib_str2=self.exp_lib_str,
                        pseudocount=self.pseudocount,
                    )
                else:
                    (
                        test_obs,
                        mean1,
                        mean2,
                        log2FC,
                        pval_ltail,
                        pval_utail,
                        pval_2tail,
                        testlist,
                    ) = stat_tools.resampling(
                        data1,
                        data2,
                        S=self.samples,
                        testFunc=stat_tools.F_mean_diff_flat,
                        permFunc=stat_tools.F_shuffle_flat,
                        adaptive=self.adaptive,
                        lib_str1=self.ctrl_lib_str,
                        lib_str2=self.exp_lib_str,
                        pseudocount=self.pseudocount,
                    )

            if self.do_histogram:
                import matplotlib.pyplot as plt

                if testlist:
                    n, bins, patches = plt.hist(
                        testlist, density=1, facecolor="c", alpha=0.75, bins=100
                    )
                else:
                    n, bins, patches = plt.hist(
                        [0, 0], density=1, facecolor="c", alpha=0.75, bins=100
                    )
                plt.xlabel("Delta Mean")
                plt.ylabel("Probability")
                plt.title("%s - Histogram of Delta Mean" % gene.orf)
                plt.axvline(test_obs, color="r", linestyle="dashed", linewidth=3)
                plt.grid(True)
                genePath = os.path.join(histPath, gene.orf + ".png")
                if not os.path.exists(histPath):
                    os.makedirs(histPath)
                plt.savefig(genePath)
                plt.clf()

            sum1 = numpy.sum(data1)
            sum2 = numpy.sum(data2)
            data.append(
                [
                    gene.orf,
                    gene.name,
                    gene.desc,
                    gene.n,
                    mean1,
                    mean2,
                    sum1,
                    sum2,
                    test_obs,
                    log2FC,
                    pval_2tail,
                ]
            )

            # Update progress
            text = "Running Resampling Method... %5.1f%%" % (100.0 * count / N)
            self.progress_update(text, count)

        #
        transit_tools.log("")  # Printing empty line to flush stdout
        transit_tools.log("Performing Benjamini-Hochberg Correction")
        data.sort()
        qval = stat_tools.BH_fdr_correction([row[-1] for row in data])

        return (data, qval)


@transit_tools.ResultsFile
class File(Analysis):
    @staticmethod
    def can_load(path):
        with open(path) as in_file:
            for line in in_file:
                if line.startswith("#"):
                    if line.startswith(Analysis.identifier):
                        return True
                else:
                    return False
        return False
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Analysis.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title="Anova",heading="",column_names=self.column_names,rows=self.rows).Show(),
                "Display Heatmap": lambda *args: self.create_heatmap(infile=self.path, output_path=self.path+".heatmap.png"),
            })
        )
        
        # 
        # get column names
        # 
        comments, headers, rows = csv.read(self.path, seperator="\t", skip_empty_lines=True)
        transit_tools.log(f'''comments = {comments}''')
        if len(comments) == 0:
            raise Exception(f'''No comments in file, and I expected the last comment to be the column names, while to load Anova file "{self.path}"''')
        self.column_names = comments[-1].split("\t")
        transit_tools.log(f'''self.column_names = {self.column_names}''')
        
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
    
    def display_histogram(self, display_frame, event):
        pass
        # gene = display_frame.grid.GetCellValue(display_frame.row, 0)
        # filepath = os.path.join(
        #     ntpath.dirname(display_frame.path),
        #     transit_tools.fetch_name(display_frame.path),
        # )
        # filename = os.path.join(filepath, gene + ".png")
        # if os.path.exists(filename):
        #     imgWindow = pytransit.file_display.ImgFrame(None, filename)
        #     imgWindow.Show()
        # else:
        #     transit_tools.show_error_dialog("Error Displaying File. Histogram image not found. Make sure results were obtained with the histogram option turned on.")
        #     print("Error Displaying File. Histogram image does not exist.")

    def create_heatmap(self, infile, output_path, topk=-1, qval=0.05, low_mean_filter=5):
        if not HAS_R:
            raise Exception(f'''Error: R and rpy2 (~= 3.0) required to run Heatmap''')
        headers = None
        data, hits = [], []
        number_of_conditions = -1

        for line in open(infile):
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


class ResamplingFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#Resampling", columns)

    def getHeader(self, path):
        DE = 0
        poslogfc = 0
        neglogfc = 0
        for line in open(path):
            if line.startswith("#"):
                continue
            tmp = line.strip().split("\t")
            if float(tmp[-1]) < 0.05:
                DE += 1
                if float(tmp[-3]) > 0:
                    poslogfc += 1
                else:
                    neglogfc += 1

        text = """Results:
    Conditionally - Essentials: %s
        Less Essential in Experimental datasets: %s
        More Essential in Experimental datasets: %s
            """ % (
            DE,
            poslogfc,
            neglogfc,
        )
        return text

    def getMenus(self):
        menus = []
        menus.append(("Display in Track View", self.displayInTrackView))
        menus.append(("Display Histogram", self.displayHistogram))
        return menus

    def displayHistogram(self, displayFrame, event):
        gene = displayFrame.grid.GetCellValue(displayFrame.row, 0)
        filepath = os.path.join(
            ntpath.dirname(displayFrame.path),
            transit_tools.fetch_name(displayFrame.path),
        )
        filename = os.path.join(filepath, gene + ".png")
        if os.path.exists(filename):
            imgWindow = pytransit.file_display.ImgFrame(None, filename)
            imgWindow.Show()
        else:
            transit_tools.show_error_dialog("Error Displaying File. Histogram image not found. Make sure results were obtained with the histogram option turned on.")
            print("Error Displaying File. Histogram image does not exist.")




    
Method = GUI = Analysis
Analysis() # make sure theres one instance