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
from pytransit.methods.analysis.pathway_enrichment_gui import Analysis as PathwayEnrichment

from pytransit.basics.lazy_dict import LazyDict

from pytransit.methods import analysis_base as base
from pytransit.tools.transit_tools import wx, pub, basename, HAS_R, FloatVector, DataFrame, StrVector, EOL
from pytransit.tools.tnseq_tools import Wig
import pytransit
import pytransit.components.file_display as file_display
import pytransit.basics.csv as csv
import pytransit.components.results_area as results_area
from pytransit.tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.universal_data import universal
from pytransit.components import parameter_panel
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.panel_helpers import make_panel, create_run_button, create_normalization_input, create_reference_condition_input, create_include_condition_list_input, create_exclude_condition_list_input, create_n_terminus_input, create_c_terminus_input, create_pseudocount_input, create_winsorize_input, create_alpha_input, create_button, create_text_box_getter, create_button, create_check_box_getter, create_control_condition_input, create_experimental_condition_input, create_preview_loess_button
command_name = sys.argv[0]

class Analysis:
    identifier  = "Resampling"
    short_name = "resampling"
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
        ctrldata=None,
        expdata=None,
        annotation_path=None,
        output_path=None,
        normalization="TTR",
        samples=10000,
        adaptive=False,
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
        do_histogram=False,
    )
    
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
        from pytransit.components.panel_helpers import Panel
        with Panel() as (self.panel, main_sizer):
            self.value_getters = LazyDict()
            sample_getter          = create_text_box_getter(self.panel, main_sizer, label_text="Samples", default_value="10000", tooltip_text="Number of samples to take when estimating the resampling histogram. More samples give more accurate estimates of the p-values at the cost of computation time.")
            
            self.value_getters.ctrldata               = create_control_condition_input(self.panel, main_sizer)
            self.value_getters.expdata                = create_experimental_condition_input(self.panel, main_sizer)
            self.value_getters.samples                = lambda *args: int(sample_getter(*args))
            self.value_getters.n_terminus             = create_n_terminus_input(self.panel, main_sizer)
            self.value_getters.c_terminus             = create_c_terminus_input(self.panel, main_sizer)
            self.value_getters.pseudocount            = create_pseudocount_input(self.panel, main_sizer)
            self.value_getters.normalization          = create_normalization_input(self.panel, main_sizer)
            self.value_getters.genome_positional_bias = create_check_box_getter(self.panel, main_sizer, label_text="Correct for Genome Positional Bias", default_value=False, tooltip_text="Check to correct read-counts for possible regional biase using LOESS. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
            create_preview_loess_button(self.panel, main_sizer)
            self.value_getters.adaptive                = create_check_box_getter(self.panel, main_sizer, label_text="Adaptive Resampling (Faster)", default_value=True, tooltip_text="Dynamically stops permutations early if it is unlikely the ORF will be significant given the results so far. Improves performance, though p-value calculations for genes that are not differentially essential will be less accurate.")
            self.value_getters.do_histogram            = create_check_box_getter(self.panel, main_sizer, label_text="Generate Resampling Histograms", default_value=False, tooltip_text="Creates .png images with the resampling histogram for each of the ORFs. Histogram images are created in a folder with the same name as the output file.")
            self.value_getters.include_zeros           = create_check_box_getter(self.panel, main_sizer, label_text="Include sites with all zeros", default_value=True, tooltip_text="Includes sites that are empty (zero) across all datasets. Unchecking this may be useful for tn5 datasets, where all nucleotides are possible insertion sites and will have a large number of empty sites (significantly slowing down computation and affecting estimates).")
            
            create_run_button(self.panel, main_sizer, from_gui_function=self.from_gui)
        
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
            transit_tools.validate_annotation(Analysis.inputs.annotation_path)
            
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
                default_file_name="resampling_output.dat",
                output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
            )
            if not Analysis.inputs.output_path:
                return None
            
            # 
            # extract universal data
            # 
            cwig_path     = universal.session_data.combined_wigs[0].main_path
            metadata_path = universal.session_data.combined_wigs[0].metadata.path
            
            from pytransit.components.samples_area import sample_table
            Analysis.inputs.combined_wig_params = dict(
                combined_wig=cwig_path,
                samples_metadata=metadata_path,
                conditions=[
                    Analysis.inputs.ctrldata,
                    Analysis.inputs.expdata,
                ],
            )
            assert Analysis.inputs.ctrldata != "[None]", "Control group can't be None"
            assert Analysis.inputs.expdata != "[None]", "Experimental group can't be None"
            
            # backwards compatibility
            Analysis.inputs.ctrldata = [Analysis.inputs.combined_wig_params["conditions"][0]]
            Analysis.inputs.expdata = [Analysis.inputs.combined_wig_params["conditions"][1]]

            Analysis.inputs.update(dict(
                annotation_path_exp=Analysis.inputs.annotation_path_exp if Analysis.inputs.diff_strains else Analysis.inputs.annotation_path
            ))
            
            return Analysis.instance

    @classmethod
    def from_args(cls, args, kwargs):
        is_combined_wig = True if kwargs.get("c", False) else False
        combined_wig_params = None
        if is_combined_wig:
            if len(args) != 5:
                print("Error: Incorrect number of args. See usage")
                print(cls.usage_string)
                sys.exit(0)
            combined_wig_params = {
                "combined_wig": kwargs.get("c"),
                "samples_metadata": args[0],
                "conditions": [args[1], args[2]],
            }
            annot_paths = args[3].split(",")
            # to show contrasted conditions for combined_wigs in output header
            ctrldata = [combined_wig_params["conditions"][0]]
            expdata = [combined_wig_params["conditions"][1]]
            output_path = args[4]
        else:
            if len(args) != 4:
                print("Error: Incorrect number of args. See usage")
                print(cls.usage_string)
                sys.exit(0)
            ctrldata = args[0].split(",")
            expdata = args[1].split(",")
            annot_paths = args[2].split(",")
            output_path = args[3]
        annotation_path = annot_paths[0]
        diff_strains = False
        annotation_path_exp = ""
        if len(annot_paths) == 2:
            annotation_path_exp = annot_paths[1]
            diff_strains = True
        if diff_strains and is_combined_wig:
            print("Error: Cannot have combined wig and different annotation files.")
            sys.exit(0)
        winz = True if "winz" in kwargs else False

        # check for unrecognized flags
        console_tools.handle_unrecognized_flags(
            "-c -s -n -h -a -ez -PC -l -iN -iC --ctrl_lib --exp_lib -Z -winz".split(),
            kwargs,
            Analysis.usage_string,
        )

        normalization = kwargs.get("n", "TTR")
        samples = int(kwargs.get("s", 10000))
        adaptive = kwargs.get("a", False)
        replicates = kwargs.get("r", "Sum")
        do_histogram = kwargs.get("h", False)
        exclude_zeros = kwargs.get("ez", False)
        include_zeros = not exclude_zeros
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
        
        cls.inputs.update(dict(
            ctrldata=ctrldata,
            expdata=expdata,
            output_path=output_path,
            normalization=normalization,
            samples=samples,
            adaptive=adaptive,
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
            annotation_path=annotation_path,
            annotation_path_exp=annotation_path,
            combined_wig_params=combined_wig_params,
            do_histogram=do_histogram,
        ))
        return Analysis.instance

    def Run(self):
        with gui_tools.nice_error_log:
            if self.inputs.do_histogram:
                try: import matplotlib.pyplot as plt
                except:
                    print("Error: cannot do histograms")
                    self.inputs.do_histogram = False

            logging.log("Starting resampling Method")
            start_time = time.time()
            if self.inputs.winz:
                logging.log("Winsorizing insertion counts")

            # Get orf data
            logging.log("Getting Data")
            if self.inputs.diff_strains:
                logging.log("Multiple annotation files found")
                logging.log(
                    "Mapping ctrl data to {0}, exp data to {1}".format(
                        self.inputs.annotation_path, self.inputs.annotation_path_exp
                    )
                )
            # 
            # Combine 
            # 
            if self.inputs.combined_wig_params:
                (position, data, filenames_in_comb_wig) = tnseq_tools.read_combined_wig(
                    self.inputs.combined_wig_params["combined_wig"]
                )
                conditions_by_file, _, _, _ = tnseq_tools.read_samples_metadata(
                    self.inputs.combined_wig_params["samples_metadata"]
                )
                condition_names = self.wigs_to_conditions(conditions_by_file, filenames_in_comb_wig)
                datasets, conditions_per_dataset = self.filter_wigs_by_conditions(
                    data, condition_names, self.inputs.combined_wig_params["conditions"],
                )
                control_condition = self.inputs.combined_wig_params["conditions"][0]
                experimental_condition = self.inputs.combined_wig_params["conditions"][1]
                data_ctrl = numpy.array(
                    [
                        each_dataset
                            for i, each_dataset in enumerate(datasets)
                                if conditions_per_dataset[i] == control_condition
                    ]
                )
                data_exp = numpy.array(
                    [
                        each_dataset
                            for i, each_dataset in enumerate(datasets)
                                if conditions_per_dataset[i] == experimental_condition
                    ]
                )
                position_ctrl, position_exp = position, position
            else:
                output = transit_tools.get_validated_data(
                    self.inputs.ctrldata, wxobj=self.wxobj
                )
                (data_ctrl, position_ctrl, *_) = output
                (data_exp, position_exp) = transit_tools.get_validated_data(
                    self.inputs.expdata, wxobj=self.wxobj
                )
            (K_ctrl, N_ctrl) = data_ctrl.shape
            (K_exp, N_exp) = data_exp.shape

            if not self.inputs.diff_strains and (N_ctrl != N_exp):
                # NOTE: returned ([], [])
                logging.error("Error: Ctrl and Exp wig files don't have the same number of sites. Make sure all .wig files come from the same strain.")

            logging.log("Preprocessing Ctrl data...")
            data_ctrl = self.preprocess_data(position_ctrl, data_ctrl)

            logging.log("Preprocessing Exp data...")
            data_exp = self.preprocess_data(position_exp, data_exp)
            
            G_ctrl = tnseq_tools.Genes(
                self.inputs.ctrldata,
                self.inputs.annotation_path,
                ignore_codon=self.inputs.ignore_codon,
                n_terminus=self.inputs.n_terminus,
                c_terminus=self.inputs.c_terminus,
                data=data_ctrl,
                position=position_ctrl,
            )
            G_exp = tnseq_tools.Genes(
                self.inputs.expdata,
                self.inputs.annotation_path_exp or self.inputs.annotation_path,
                ignore_codon=self.inputs.ignore_codon,
                n_terminus=self.inputs.n_terminus,
                c_terminus=self.inputs.c_terminus,
                data=data_exp,
                position=position_exp,
            )

            do_library_resampling = False
            # If library string not empty
            if self.inputs.ctrl_lib_str or self.inputs.exp_lib_str:
                letters_ctrl = set(self.inputs.ctrl_lib_str)
                letters_exp = set(self.inputs.exp_lib_str)

                # Check if using exactly 1 letters; i.e. no different libraries
                if len(letters_ctrl) == 1 and letters_exp == 1:
                    pass
                # If using more than one letter, then check no differences in set
                else:
                    lib_diff = letters_ctrl ^ letters_exp
                    # Check that their differences
                    if not lib_diff:
                        do_library_resampling = True
                    else:
                        # NOTE: kept going (set self.inputs.ctrl_lib_str = "", self.inputs.exp_lib_str = "")
                        logging.error(
                            "Error: Library Strings (Ctrl = %s, Exp = %s) do not use the same letters. Make sure every letter / library is represented in both Control and Experimental Conditions. Proceeding with resampling assuming all datasets belong to the same library."
                            % (self.inputs.ctrl_lib_str, self.inputs.exp_lib_str)
                        )
            
            (data, qval) = self.run_resampling(G_ctrl, G_exp, do_library_resampling)
            # 
            # write output
            # 
            if True:
                # 
                # generate rows
                # 
                rows = []
                for row_index, row in enumerate(data):
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
                    if self.inputs.Z == True:
                        p = pval_2tail / 2  # convert from 2-sided back to 1-sided
                        if p == 0:
                            p = 1e-5  # or 1 level deeper the num of iterations of resampling, which is 1e-4=1/10000, by default
                        if p == 1:
                            p = 1 - 1e-5
                        z = scipy.stats.norm.ppf(p)
                        if log2FC > 0:
                            z *= -1
                        rows.append(
                            (
                                "%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.1f\t%1.2f\t%1.1f\t%1.5f\t%0.2f\t%1.5f"
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
                                    qval[row_index],
                                )
                            ).split('\t')
                        )
                    else:
                        rows.append(
                            (
                                "%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.1f\t%1.2f\t%1.1f\t%1.5f\t%1.5f"
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
                                    qval[row_index],
                                )
                            ).split('\t')
                        )
                
                # 
                # write to file
                # 
                transit_tools.write_result(
                    path=self.inputs.output_path,
                    file_kind=Analysis.identifier,
                    rows=rows,
                    column_names=Analysis.columns if not self.inputs.Z else [
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
                    ],
                    extra_info=dict(
                        parameters=dict(
                            samples=self.inputs.samples,
                            norm=self.inputs.normalization,
                            histograms=self.inputs.do_histogram,
                            adaptive=self.inputs.adaptive,
                            exclude_zeros=not self.inputs.include_zeros,
                            pseudocounts=self.inputs.pseudocount,
                            LOESS=self.inputs.LOESS,
                            n_terminus=self.inputs.n_terminus,
                            c_terminus=self.inputs.c_terminus,
                        ),
                        control_data=(",".join(self.inputs.ctrldata)),
                        experimental_data=(",".join(self.inputs.expdata)),
                        annotation_path=self.inputs.annotation_path,
                        **({} if not self.inputs.diff_strains else dict(
                            annotation_path_exp=self.inputs.annotation_path_exp,
                        )),
                        time=(time.time() - start_time),
                    ),
                )
                results_area.add(self.inputs.output_path)
                
            logging.log(f"Finished running {Analysis.short_name}")

    def preprocess_data(self, position, data):
        (K, N) = data.shape

        if self.inputs.normalization != "nonorm":
            logging.log("Normalizing using: %s" % self.inputs.normalization)
            (data, factors) = norm_tools.normalize_data(
                data,
                self.inputs.normalization,
                self.inputs.ctrldata + self.inputs.expdata,
                self.inputs.annotation_path,
            )

        if self.inputs.LOESS:
            logging.log("Performing LOESS Correction")
            from pytransit.tools import stat_tools
            for j in range(K):
                data[j] = stat_tools.loess_correction(position, data[j])

        return data

    def wigs_to_conditions(self, conditions_by_file, filenames_in_comb_wig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditions_by_file.get(f, None) for f in filenames_in_comb_wig]

    def filter_wigs_by_conditions(self, data, conditions, included_conditions):
        """
            Filters conditions from wig to ctrl, exp conditions only
            ([[Wigdata]], [ConditionCtrl, ConditionExp]) -> Tuple([[Wigdata]], [Condition])
        """
        d_filtered, cond_filtered = [], []
        if len(included_conditions) != 2:
            logging.error("Only 2 conditions expected, but got:", included_conditions)
        
        for i, c in enumerate(conditions):
            if c in included_conditions:
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])

        return (numpy.array(d_filtered), numpy.array(cond_filtered))

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
        self, G_ctrl, G_exp=None, do_library_resampling=False
    ):
        from pytransit.tools import stat_tools
        
        data = []
        N = len(G_ctrl)
        count = 0
        
        if self.inputs.do_histogram:
            import matplotlib.pyplot as plt
            hist_path = os.path.join(
                os.path.dirname(self.inputs.output_path),
                transit_tools.fetch_name(self.inputs.output_path) + "_histograms",
            )
            os.makedirs(hist_path, exist_ok=True)
        
        for gene in G_ctrl:
            if gene.orf not in G_exp:
                if self.inputs.diff_strains:
                    continue
                else:
                    # NOTE: returned ([], [])
                    logging.error("Error: Gene in ctrl data not present in exp data. Make sure all .wig files come from the same strain.")

            gene_exp = G_exp[gene.orf]
            count += 1
            
            if not self.inputs.diff_strains and gene.n != gene_exp.n:
                # NOTE: returned ([], [])
                logging.error("Error: No. of TA sites in Exp and Ctrl data are different. Make sure all .wig files come from the same strain.")

            if (gene.k == 0 and gene_exp.k == 0) or gene.n == 0 or gene_exp.n == 0:
                test_obs   = 0
                mean1      = 0
                mean2      = 0
                log2FC     = 0
                pval_ltail = 1.00
                pval_utail = 1.00
                pval_2tail = 1.00
                testlist   = []
                data1      = [0]
                data2      = [0]
            else:
                if not self.inputs.include_zeros:
                    ii_ctrl = numpy.sum(gene.reads, axis=0) > 0
                    ii_exp = numpy.sum(gene_exp.reads, axis=0) > 0
                else:
                    ii_ctrl = numpy.ones(gene.n) == 1
                    ii_exp = numpy.ones(gene_exp.n) == 1

                # data1 = gene.reads[:,ii_ctrl].flatten() + self.inputs.pseudocount # we used to have an option to add pseudocounts to each observation, like this
                data1 = gene.reads[:, ii_ctrl].flatten()
                data2 = gene_exp.reads[:, ii_exp].flatten()
                if self.inputs.winz:
                    data1 = self.winsorize_resampling(data1)
                    data2 = self.winsorize_resampling(data2)
                
                if do_library_resampling:
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
                        S=self.inputs.samples,
                        test_func=stat_tools.f_mean_diff_dict,
                        perm_func=stat_tools.f_shuffle_dict_libraries,
                        adaptive=self.inputs.adaptive,
                        lib_str1=self.inputs.ctrl_lib_str,
                        lib_str2=self.inputs.exp_lib_str,
                        pseudocount=self.inputs.pseudocount,
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
                        S=self.inputs.samples,
                        test_func=stat_tools.f_mean_diff_flat,
                        perm_func=stat_tools.f_shuffle_flat,
                        adaptive=self.inputs.adaptive,
                        lib_str1=self.inputs.ctrl_lib_str,
                        lib_str2=self.inputs.exp_lib_str,
                        pseudocount=self.inputs.pseudocount,
                    )
                
                # 
                # write historgram if needed
                # 
                if self.inputs.do_histogram:
                    # TODO: this should probably be done at the bottom of .Run() instead of inside .run_resampling()
                    import matplotlib.pyplot as plt
                    
                    testlist = testlist or [0,0]
                    n, bins, patches = plt.hist(testlist, density=1, facecolor="c", alpha=0.75, bins=100)
                    plt.xlabel("Delta Mean")
                    plt.ylabel("Probability")
                    plt.title("%s - Histogram of Delta Mean" % gene.orf)
                    plt.axvline(test_obs, color="r", linestyle="dashed", linewidth=3)
                    plt.grid(True)
                    gene_path = os.path.join(hist_path, gene.orf + ".png")
                    plt.savefig(gene_path)
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
            percentage = (100.0 * count / N)
            text = "Running Resampling Method... %5.1f%%" % percentage
            parameter_panel.progress_update(text, percentage)

        logging.log("")  # Printing empty line to flush stdout
        logging.log("Performing Benjamini-Hochberg Correction")
        data.sort()
        qval = stat_tools.bh_fdr_correction([row[-1] for row in data])

        return (data, qval)


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
                "Display Volcano Plot": lambda *args: self.graph_volcano_plot(),
                "Pathway Enrichment": lambda *args: PathwayEnrichment.call_from_results_panel(path),
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
    
    def graph_volcano_plot(self):
        # questions:
            # are the selected rows correct ("log2FC", "Adj. p-value")?
            # what is the q_value supposed to be?
            # why are some log2 and the other axis log10?
        with gui_tools.nice_error_log:
            try: import matplotlib.pyplot as plt
            except:
                print("Error: cannot do plots, no matplotlib")
                
            log2_fc_values = [ each_row["log2FC"]  for each_row in self.rows ]
            p_values       = [ each_row["p-value"] for each_row in self.rows ]
            q_values       = [ each_row["Adj. p-value"] for each_row in self.rows ]
            log10_p_values = []
            for each_p_value in p_values:
                try:
                    log10_p_value = -math.log(float(each_p_value), 10)
                except ValueError as e:
                    log10_p_value = None
                
                log10_p_values.append(log10_p_value)
            
            # 
            # replace missing values
            # 
            good_values = [ each for each in log10_p_values if each is not None ]
            max_log10_p_values = max(good_values)
            # replace None values with max value
            log10_p_values = [ (each if each is not None else max_log10_p_values) for each in log10_p_values ]
            
            # 
            # compute threshold
            # 
            # NOTE: find the p-value (horizontal line) that corrisponds to where q-value == 0.05
            if True:
                threshold = 0.00001
                backup_thresh = 0.00001
                q_and_p_values_list = list(zip(q_values, p_values))
                q_and_p_values_list.sort()
                for (q, p) in q_and_p_values_list:
                    backup_thresh = p
                    if q > 0.05:
                        break
                    threshold = p # should be equivlent to: threshold = max(p, threshold) # because they're sorted

                if threshold == 0:
                    threshold = backup_thresh
            
            # 
            # plot (log2_fc_values, log10_p_values, threshold)
            # 
            plt.plot(log2_fc_values, log10_p_values, "bo")
            plt.axhline( -math.log(threshold, 10), color="r", linestyle="dashed", linewidth=3)
            plt.xlabel("Log Fold Change (base 2)")
            plt.ylabel("-Log p-value (base 10)")
            plt.suptitle("Resampling - Volcano plot")
            plt.title("Adjusted threshold (red line): P-value=%1.8f" % threshold)
            plt.show()

class ResamplingFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#Resampling", columns)

    def get_header(self, path):
        DE = 0
        poslogfc = 0
        neglogfc = 0
        with open(path) as file:
            for line in file:
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

    def get_menus(self):
        menus = []
        menus.append(("Display in Track View", self.display_in_track_view))
        return menus

    
Method = GUI = Analysis
Analysis() # make sure there's one instance
