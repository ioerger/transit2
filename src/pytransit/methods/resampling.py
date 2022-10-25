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
from pytransit.components.parameter_panel import set_instructions
from pytransit.methods.pathway_enrichment import Method as PathwayEnrichment

from pytransit.generic_tools.lazy_dict import LazyDict

from pytransit.specific_tools.transit_tools import wx, basename, HAS_R, FloatVector, DataFrame, StrVector
from pytransit.specific_tools.tnseq_tools import Wig
import pytransit
import pytransit.components.file_display as file_display
from pytransit.generic_tools import csv, misc, informative_iterator
import pytransit.components.results_area as results_area
from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components import parameter_panel
from pytransit.components.spreadsheet import SpreadSheet


@misc.singleton
class Method:
    name        = "Resampling"
    identifier  = name
    cli_name    = identifier.lower()
    menu_name   = f"{identifier} - test of conditional essentiality between two conditions"
    description = """Method for determining conditional essentiality based on resampling (i.e. permutation test). Identifies significant changes in mean read-counts for each gene after normalization."""
    significance_threshold = 0.05
    
    column_names = [
        "ORF",
        "Gene Name",
        "Description",
        "Sites",
        "Mean Control",
        "Mean Experimental",
        "Log 2 FC",
        "Sum Control",
        "Sum Experimental",
        "Delta Mean",
        "P Value",
        "Adj P Value",
    ]
    
    valid_cli_flags = [
        "-c",
        "-s",
        "-n",
        "-h",
        "-a",
        "-ez",
        "-PC",
        "-l",
        "-iN",
        "-iC",
        "--ctrl_lib",
        "--exp_lib",
        "-Z",
        "-winz",
        "-sr",
    ]
    
    inputs = LazyDict(
        ctrldata=None,
        expdata=None,
        annotation_path=None,
        output_path=None,
        normalization="TTR",
        samples=10000,
        adaptive=False,
        include_zeros=True, # FIXME: include_zeros=False breaks resampling (even on master I believe)
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
        site_restricted=False,
    )
    
    usage_string = f"""
        {console_tools.subcommand_prefix} resampling <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        ---
        OR
        ---
        {console_tools.subcommand_prefix} resampling -c <combined wig file> <samples_metadata file> <ctrl condition name> <exp condition name> <annotation .prot_table> <output file> [Optional Arguments]
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
        -sr             :=  site-restricted resampling; more sensitive, might find a few more significant conditionally essential genes"
    """.replace("\n        ", "\n")
    
    @gui.add_menu("Method", "himar1", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    @gui.add_menu("Method", "tn5", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)

    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            set_instructions(
                method_short_text=self.name,
                method_long_text="",
                method_specific_instructions="""
                The resampling method is a comparative analysis the allows that can be used to determine conditional essentiality of genes. It is based on a permutation test, and is capable of determining read-counts that are significantly different across conditions.

                See Pathway Enrichment Method for post-processing the hits to determine if the hits are associated with a particular functional catogory of genes or known biological pathway.
                
                1. Of the Conditions in the Conditions pane, select one to be the control condition using the 'Control Condition' dropdown

                2. Of the Conditions in the Conditions pane, select one to be the experimental condition using the 'Experimental Condition' dropdown

                3.[Optional] Select/Adjust other parameters

                4.[Optional] Select from the samples panel and then click on 'Preview LOESS fit' to see the loess fit graph. This is the equivalent of selecting values from the samples panel and selecting 'LOESS' on the dropdown

                5.[Optional] If you select to 'Generate Resampling Histograms', a folder titled 'resampling_output_histograms' will be generated and populated locally

                6. Click Run
                """.replace("\n                    ","\n"),
            )

            self.value_getters = LazyDict()
            
            self.value_getters.ctrldata               = panel_helpers.create_control_condition_input(panel, main_sizer)
            self.value_getters.expdata                = panel_helpers.create_experimental_condition_input(panel, main_sizer)
            self.value_getters.samples                = panel_helpers.create_int_getter(panel, main_sizer, label_text="Samples", default_value="10000", tooltip_text="Number of samples to take when estimating the resampling histogram. More samples give more accurate estimates of the p-values at the cost of computation time.")
            self.value_getters.n_terminus             = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus             = panel_helpers.create_c_terminus_input(panel, main_sizer)
            self.value_getters.pseudocount            = panel_helpers.create_pseudocount_input(panel, main_sizer)
            self.value_getters.normalization          = panel_helpers.create_normalization_input(panel, main_sizer)
            self.value_getters.site_restricted        = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="Site-restricted resampling", default_value=False, tooltip_text="Restrict permutations of insertion counts in a gene to each individual TA site, which could be more sensitive (detect more conditional-essentials) than permuting counts over all TA sites pooled (which is the default).")
            self.value_getters.genome_positional_bias = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="Correct for Genome Positional Bias", default_value=False, tooltip_text="Check to correct read-counts for possible regional biase using LOESS. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
            panel_helpers.create_preview_loess_button(panel, main_sizer) # TODO: change tooltip text that references this button, then remove this button
            self.value_getters.adaptive                = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="Adaptive Resampling (Faster)", default_value=True, tooltip_text="Dynamically stops permutations early if it is unlikely the ORF will be significant given the results so far. Improves performance, though p-value calculations for genes that are not differentially essential will be less accurate.")
            self.value_getters.do_histogram            = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="Generate Resampling Histograms", default_value=False, tooltip_text="Creates .png images with the resampling histogram for each of the ORFs. Histogram images are created in a folder with the same name as the output file.")
            self.value_getters.include_zeros           = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="Include sites with all zeros", default_value=True, tooltip_text="Includes sites that are empty (zero) across all datasets. Unchecking this may be useful for tn5 datasets, where all nucleotides are possible insertion sites and will have a large number of empty sites (significantly slowing down computation and affecting estimates).")
            
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
        
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
        # save result files
        # 
        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}_output.txt",
            output_extensions='Common output extensions (*.csv,*.dat,*.txt,*.out)|*.csv;*.dat;*.txt;*.out;|\nAll files (*.*)|*.*',
        )
        if not Method.inputs.output_path:
            return None
        
        # 
        # extract universal data
        # 
        cwig_path     = gui.combined_wigs[-1].main_path
        metadata_path = gui.combined_wigs[-1].metadata.path
        
        Method.inputs.combined_wig_params = dict(
            combined_wig=cwig_path,
            samples_metadata=metadata_path,
            conditions=[
                Method.inputs.ctrldata,
                Method.inputs.expdata,
            ],
        )
        assert Method.inputs.ctrldata != "[None]", "Control group can't be None"
        assert Method.inputs.expdata != "[None]", "Experimental group can't be None"
        
        # backwards compatibility
        Method.inputs.ctrldata = [Method.inputs.combined_wig_params["conditions"][0]]
        Method.inputs.expdata = [Method.inputs.combined_wig_params["conditions"][1]]

        Method.inputs.update(dict(
            annotation_path_exp=Method.inputs.annotation_path_exp if Method.inputs.diff_strains else Method.inputs.annotation_path
        ))
        
        return Method

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        
        is_combined_wig = True if kwargs.get("c", False) else False
        combined_wig_params = None
        if is_combined_wig:
            console_tools.enforce_number_of_args(args, Method.usage_string, exactly=5)
            combined_wig_params = {
                "combined_wig": kwargs.get("c"),
                "samples_metadata": args[0],
                "conditions": [args[1], args[2]],
            }
            annot_paths = args[3].split(",")
            # to show contrasted conditions for combined_wigs in output header
            ctrldata = [combined_wig_params["conditions"][0]]
            expdata  = [combined_wig_params["conditions"][1]]
            output_path = args[4]
        else:
            console_tools.enforce_number_of_args(args, Method.usage_string, exactly=4)
            ctrldata    = args[0].split(",")
            expdata     = args[1].split(",")
            annot_paths = args[2].split(",")
            output_path = args[3]
        
        annotation_path = annot_paths[0]
        diff_strains = False
        annotation_path_exp = ""
        if len(annot_paths) == 2:
            annotation_path_exp = annot_paths[1]
            diff_strains = True
        if diff_strains and is_combined_wig:
            logging.error("Error: Cannot have combined wig and different annotation files.")

        winz          = True if "winz" in kwargs else False
        normalization = kwargs.get("n", Method.inputs.normalization)
        samples       = int(kwargs.get("s", Method.inputs.samples))
        adaptive      = kwargs.get("a", Method.inputs.adaptive)
        replicates    = kwargs.get("r", Method.inputs.replicates)
        do_histogram  = kwargs.get("h", Method.inputs.do_histogram)
        include_zeros = not kwargs.get("ez", not Method.inputs.include_zeros)
        pseudocount   = float(kwargs.get("PC", Method.inputs.pseudocount))  # use -PC (new semantics: for LFCs) instead of -pc (old semantics: fake counts)
        
        Z = True if "Z" in kwargs else False

        LOESS = kwargs.get("l", False)
        ignore_codon = True

        n_terminus = float(kwargs.get("iN", 0.00))  # integer interpreted as percentage
        c_terminus = float(kwargs.get("iC", 0.00))
        ctrl_lib_str = kwargs.get("-ctrl_lib", "")
        exp_lib_str = kwargs.get("-exp_lib", "")
        
        Method.inputs.update(dict(
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
        Method.Run()

    def Run(self):
        if self.inputs.do_histogram:
            try: import matplotlib.pyplot as plt
            except:
                print("Error: cannot do histograms")
                self.inputs.do_histogram = False

        logging.log("Starting resampling Method")
        start_time = time.time()
        if self.inputs.winz:
            logging.log("Winsorizing insertion counts")
        
        if self.inputs.ctrl_lib_str and self.inputs.site_restricted:
            raise Exception("Cannot do site_restricted resampling with library strings at same time")

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
            (position, data, filenames_in_comb_wig) = tnseq_tools.CombinedWigData.load(
                self.inputs.combined_wig_params["combined_wig"]
            )
            conditions_by_wig_fingerprint, _, _, _ = tnseq_tools.read_samples_metadata(
                self.inputs.combined_wig_params["samples_metadata"]
            )
            condition_names = self.wigs_to_conditions(conditions_by_wig_fingerprint, filenames_in_comb_wig)
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
            output = transit_tools.get_validated_data(self.inputs.ctrldata)
            (data_ctrl, position_ctrl, *_) = output
            (data_exp, position_exp) = transit_tools.get_validated_data(self.inputs.expdata)
        (K_ctrl, N_ctrl) = data_ctrl.shape
        (K_exp, N_exp) = data_exp.shape

        if not self.inputs.diff_strains and (N_ctrl != N_exp):
            # NOTE: returned ([], [])
            logging.error("Error: Ctrl and Exp wig files don't have the same number of sites. Make sure all .wig files come from the same strain.")

        logging.log("Preprocessing Ctrl data...")
        data_ctrl = self.preprocess_data(position_ctrl, data_ctrl)

        logging.log("Preprocessing Exp data...")
        data_exp = self.preprocess_data(position_exp, data_exp)
        
        g_ctrl = tnseq_tools.Genes(
            self.inputs.ctrldata,
            self.inputs.annotation_path,
            ignore_codon=self.inputs.ignore_codon,
            n_terminus=self.inputs.n_terminus,
            c_terminus=self.inputs.c_terminus,
            data=data_ctrl,
            position=position_ctrl,
        )
        g_exp = tnseq_tools.Genes(
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
        
        (data, qval) = self.run_resampling(g_ctrl, g_exp, do_library_resampling)
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
                file_kind=Method.identifier,
                rows=rows,
                column_names=Method.column_names if not self.inputs.Z else [
                    "ORF",
                    "Gene Name",
                    "Description",
                    "Sites",
                    "Mean Control",
                    "Mean Experimental",
                    "Log 2 FC",
                    "Sum Control",
                    "Sum Experimental",
                    "Delta Mean",
                    "P Value",
                    "Z Score",
                    "Adj P Value",
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
                        site_restricted=self.inputs.site_restricted,
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
            
        logging.log(f"Finished running {Method.identifier}")

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
            from pytransit.specific_tools import stat_tools
            for j in range(K):
                data[j] = stat_tools.loess_correction(position, data[j])

        return data

    def wigs_to_conditions(self, conditions_by_wig_fingerprint, filenames_in_comb_wig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditions_by_wig_fingerprint.get(f, None) for f in filenames_in_comb_wig]

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

    def winsorize_for_resampling(self, data):
        """
        Arguments:
            data:
                input is a 2D array of insertion counts for gene (not pre-flattened)
        """
        
        original_shape = data.shape
        assert len(original_shape)==2, "winsorize_resampling() expected 2D numpy array"
        counts = data.flatten().tolist()
        if len(counts) < 3:
            return data
        
        s = sorted(counts, reverse=True)
        if s[1] == 0:
            return data # don't do anything if there is only 1 non-zero value
        
        c2 = [
            (s[1] if x==s[0] else x)
                for x in counts 
        ]
        
        return numpy.array(c2).reshape(original_shape)

    def run_resampling(
        self, g_ctrl, g_exp=None, do_library_resampling=False
    ):
        from pytransit.specific_tools import stat_tools
        
        data = []
        control_group_size = len(g_ctrl)
        count = 0
        
        if self.inputs.do_histogram:
            import matplotlib.pyplot as plt
            hist_path = os.path.join(
                os.path.dirname(self.inputs.output_path),
                transit_tools.fetch_name(self.inputs.output_path) + "_histograms",
            )
            os.makedirs(hist_path, exist_ok=True)
        
        for progress, gene in informative_iterator.ProgressBar(g_ctrl, title="Running Resampling "):
            if gene.orf not in g_exp:
                if self.inputs.diff_strains:
                    continue
                else:
                    # NOTE: returned ([], [])
                    logging.error("Error: Gene in ctrl data not present in exp data. Make sure all .wig files come from the same strain.")

            gene_exp = g_exp[gene.orf]
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
                data1 = gene.reads[:,ii_ctrl]###.flatten() #TRI - do not flatten, as of 9/6/22
                data2 = gene_exp.reads[:,ii_exp]###.flatten()
                if self.inputs.winz:
                    data1 = self.winsorize_for_resampling(data1)
                    data2 = self.winsorize_for_resampling(data2)
                
                #data1 = gene.reads[:,ii_ctrl].flatten() + self.pseudocount # we used to have an option to add pseudocounts to each observation, like this
                data1 = gene.reads[:,ii_ctrl]###.flatten() #TRI - do not flatten, as of 9/6/22
                data2 = gene_exp.reads[:,ii_exp]###.flatten()

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
                        site_restricted=self.inputs.site_restricted,
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
                        site_restricted=self.inputs.site_restricted,
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
            percentage = (100.0 * count / control_group_size)
            if gui.is_active:
                text = "Running Resampling Method... %5.1f%%" % percentage
                parameter_panel.progress_update(text, percentage)

        logging.log("")  # Printing empty line to flush stdout
        logging.log("Performing Benjamini-Hochberg Correction")
        data.sort()
        qval = stat_tools.bh_fdr_correction([row[-1] for row in data])

        return (data, qval)


@transit_tools.ResultsFile
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier,
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[ "Adj P Value", "P Value" ]
                ).Show(),
                "Display Volcano Plot": lambda *args: self.graph_volcano_plot(),
                "Pathway Enrichment": lambda *args: PathwayEnrichment.call_from_results_panel(path),
            })
        )
        
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        parameters = LazyDict(self.extra_data.get("parameters", {}))
        number_of_significant = len([ 1 for each_row in self.rows if each_row["Adj P Value"] < Method.significance_threshold ])
        self.values_for_result_table.update({
            " ": f"{number_of_significant} significant conditionally essential genes"
        })
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
    
    def graph_volcano_plot(self):
        # questions:
            # are the selected rows correct ("log2FC", "Adj P Value")?
            # what is the q_value supposed to be?
            # why are some log2 and the other axis log10?
            with gui_tools.nice_error_log:
                Method.inputs.volcano_output_path = gui_tools.ask_for_output_file_path(
                default_file_name=f"volcano.png",
                output_extensions='PNG file (*.png)|*.png;|\nAll files (*.*)|*.*',
            )
    
            try: import matplotlib.pyplot as plt
            except:
                print("Error: cannot do plots, no matplotlib")
                
            log2_fc_values = [ each_row["Log 2 FC"]  for each_row in self.rows ]
            p_values       = [ each_row["P Value"] for each_row in self.rows ]
            q_values       = [ each_row["Adj P Value"] for each_row in self.rows ]
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
            plt.scatter(log2_fc_values, log10_p_values, marker=".")
            plt.axhline( -math.log(threshold, 10), color="r", linestyle="dashed", linewidth=3)
            plt.xlabel("Log Fold Change (base 2)")
            plt.ylabel("-Log P Value (base 10)")
            plt.suptitle("Resampling - Volcano plot")
            plt.title("Adjusted Threshold (red line): P Value=%1.8f" % threshold)
            #plt.show()
            plt.savefig(Method.inputs.volcano_output_path, bbox_inches='tight')

            if gui.is_active:
                logging.log(f"Adding File: {Method.inputs.volcano_output_path}")
                results_area.add(Method.inputs.volcano_output_path)

