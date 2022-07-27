import sys
import os
import time
import ntpath
import math
import random
import datetime
import heapq

import numpy
import scipy.stats
from super_map import LazyDict

from pytransit.transit_tools import wx, pub
from pytransit.analysis import base
import pytransit
import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools
from pytransit.core_data import universal

default_padding = 5 # not sure what the units are
first_arg = sys.argv[0]


############# GUI ELEMENTS ##################

short_name = "anova -- test"
long_name = "AnovaGUI"
short_desc = "AnovaGUI"
long_desc = """Anova GUI"""
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
usage_string = f"""
    python3 {first_arg} test1 <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    ---
    OR
    ---
    python3 {first_arg} test1 -c <combined wig file> <samples_metadata file> <ctrl condition name> <exp condition name> <annotation .prot_table> <output file> [Optional Arguments]
    NB: The ctrl and exp condition names should match Condition names in samples_metadata file.

    Optional Arguments:
    -s <integer>    :=  Number of samples. Default: -s 10000
    -n <string>     :=  Normalization method. Default: -n TTR
    -h              :=  Output histogram of the permutations for each gene. Default: Turned Off.
    -a              :=  Perform adaptive test1. Default: Turned Off.
    -ez             :=  Exclude rows with zero across conditions. Default: Turned off
                        (i.e. include rows with zeros).
    -PC <float>     :=  Pseudocounts used in calculating LFC. (default: 1)
    -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias.
                        Default: Turned Off.
    -iN <int>       :=  Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 0
    -iC <int>       :=  Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 0
    --ctrl_lib      :=  String of letters representing library of control files in order
                        e.g. 'AABB'. Default empty. Letters used must also be used in --exp_lib
                        If non-empty, test1 will limit permutations to within-libraries.

    --exp_lib       :=  String of letters representing library of experimental files in order
                        e.g. 'ABAB'. Default empty. Letters used must also be used in --ctrl_lib
                        If non-empty, test1 will limit permutations to within-libraries.
    -winz           :=  winsorize insertion counts for each gene in each condition 
                        (replace max cnt in each gene with 2nd highest; helps mitigate effect of outliers)
""".replace("\n    ", "    ")

parameters = LazyDict(
    combined_wig        = None
    metadata            = None
    annotation          = None
    normalization       = None
    output_file         = None
    excluded_conditions = []
    included_conditions = []
    n_terminus               = 0.0
    c_terminus               = 0.0
    pseudocount         = 1
    winz                = False
    refs                = []
    alpha               = 1000
)

def when_data_has_been_collected():
    try:
        import matplotlib.pyplot as plt
    except:
        print("Error: cannot do histograms")
        self.doHistogram = False

    gui_tools.set_status("Starting test1 Method")
    start_time = time.time()
    if self.winz:
        gui_tools.set_status("Winsorizing insertion counts")

    histPath = ""
    if self.doHistogram:
        histPath = os.path.join(
            os.path.dirname(self.output.name),
            transit_tools.fetch_name(self.output.name) + "_histograms",
        )
        if not os.path.isdir(histPath):
            os.makedirs(histPath)

    # Get orf data
    gui_tools.set_status("Getting Data")
    if self.diffStrains:
        gui_tools.set_status("Multiple annotation files found")
        gui_tools.set_status(
            "Mapping ctrl data to {0}, exp data to {1}".format(
                self.annotation_path, self.annotation_path_exp
            )
        )

    if self.combinedWigParams:
        (position, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(
            self.combinedWigParams["combined_wig"]
        )
        conditionsByFile, _, _, _ = tnseq_tools.read_samples_metadata(
            self.combinedWigParams["samples_metadata"]
        )
        conditions = self.wigs_to_conditions(conditionsByFile, filenamesInCombWig)
        data, conditions = self.filter_wigs_by_conditions(
            data, conditions, self.combinedWigParams["conditions"]
        )
        data_ctrl = numpy.array(
            [
                d
                for i, d in enumerate(data)
                if conditions[i].lower() == self.combinedWigParams["conditions"][0]
            ]
        )
        data_exp = numpy.array(
            [
                d
                for i, d in enumerate(data)
                if conditions[i].lower() == self.combinedWigParams["conditions"][1]
            ]
        )
        position_ctrl, position_exp = position, position
    else:
        (data_ctrl, position_ctrl) = transit_tools.get_validated_data(
            self.ctrldata, frame=universal.frame
        )
        (data_exp, position_exp) = transit_tools.get_validated_data(
            self.expdata, frame=universal.frame
        )
    (K_ctrl, N_ctrl) = data_ctrl.shape
    (K_exp, N_exp) = data_exp.shape

    if not self.diffStrains and (N_ctrl != N_exp):
        self.transit_error(
            "Error: Ctrl and Exp wig files don't have the same number of sites."
        )
        self.transit_error("Make sure all .wig files come from the same strain.")
        return
    # (data, position) = transit_tools.get_validated_data(self.ctrldata+self.expdata, frame=universal.frame)

    gui_tools.set_status("Preprocessing Ctrl data...")
    data_ctrl = self.preprocess_data(position_ctrl, data_ctrl)

    gui_tools.set_status("Preprocessing Exp data...")
    data_exp = self.preprocess_data(position_exp, data_exp)

    G_ctrl = tnseq_tools.Genes(
        self.ctrldata,
        self.annotation_path,
        ignoreCodon=self.ignoreCodon,
        n_terminus=self.n_terminus,
        c_terminus=self.c_terminus,
        data=data_ctrl,
        position=position_ctrl,
    )
    G_exp = tnseq_tools.Genes(
        self.expdata,
        self.annotation_path_exp,
        ignoreCodon=self.ignoreCodon,
        n_terminus=self.n_terminus,
        c_terminus=self.c_terminus,
        data=data_exp,
        position=position_exp,
    )

    doLibraryAnova = False
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
                doLibraryAnova = True
            else:
                transit_tools.transit_error(
                    "Error: Library Strings (Ctrl = %s, Exp = %s) do not use the same letters. Make sure every letter / library is represented in both Control and Experimental Conditions. Proceeding with test1 assuming all datasets belong to the same library."
                    % (self.ctrl_lib_str, self.exp_lib_str)
                )
                self.ctrl_lib_str = ""
                self.exp_lib_str = ""

    (data, qval) = self.run_test1(G_ctrl, G_exp, doLibraryAnova, histPath)
    self.write_output(data, qval, start_time)

    self.finish()
    gui_tools.set_status("Finished test1 Method")

def collect_data_from_args(rawargs):
    (args, kwargs) = transit_tools.clean_args(rawargs)

    isCombinedWig = True if kwargs.get("c", False) else False
    combinedWigParams = None
    if isCombinedWig:
        if len(args) != 5:
            print("Error: Incorrect number of args. See usage")
            print(self.usage_string())
            sys.exit(0)
        combinedWigParams = {
            "combined_wig": kwargs.get("c"),
            "samples_metadata": args[0],
            "conditions": [args[1].lower(), args[2].lower()],
        }
        annot_paths = args[3].split(",")
        # to show contrasted conditions for combined_wigs in output header
        ctrldata = [combinedWigParams["conditions"][0]]
        expdata = [combinedWigParams["conditions"][1]]
        output_path = args[4]
    else:
        if len(args) != 4:
            print("Error: Incorrect number of args. See usage")
            print(self.usage_string())
            sys.exit(0)
        ctrldata = args[0].split(",")
        expdata = args[1].split(",")
        annot_paths = args[2].split(",")
        output_path = args[3]
    annotation_path = annot_paths[0]
    diffStrains = False
    annotation_path_exp = ""
    if len(annot_paths) == 2:
        annotation_path_exp = annot_paths[1]
        diffStrains = True
    if diffStrains and isCombinedWig:
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
            print(ZinbMethod.usage_string())
            sys.exit(0)

    normalization = kwargs.get("n", "TTR")
    samples       = int(kwargs.get("s", 10000))
    adaptive      = kwargs.get("a", False)
    doHistogram   = kwargs.get("h", False)
    replicates    = kwargs.get("r", "Sum")
    excludeZeros  = kwargs.get("ez", False)
    includeZeros  = not excludeZeros
    pseudocount   = float(
        kwargs.get("PC", 1.0)
    )  # use -PC (new semantics: for LFCs) instead of -pc (old semantics: fake counts)

    Z = True if "Z" in kwargs else False

    LOESS = kwargs.get("l", False)
    ignoreCodon = True

    n_terminus = float(kwargs.get("iN", 0.00))  # integer interpreted as percentage
    c_terminus = float(kwargs.get("iC", 0.00))
    ctrl_lib_str = kwargs.get("-ctrl_lib", "")
    exp_lib_str = kwargs.get("-exp_lib", "")
    
    parameters.update(dict(
        ctrldata=ctrldata,
        expdata=expdata,
        annotation_path=annotation_path,
        output_file=output_file,
        normalization=normalization,
        samples=samples,
        adaptive=adaptive,
        doHistogram=doHistogram,
        includeZeros=includeZeros,
        pseudocount=pseudocount,
        replicates=replicates,
        LOESS=LOESS,
        ignoreCodon=ignoreCodon,
        n_terminus=n_terminus,
        c_terminus=c_terminus,
        ctrl_lib_str=ctrl_lib_str,
        exp_lib_str=exp_lib_str,
        winz=winz,
        Z=Z,
        diffStrains=diffStrains,
        annotation_path_exp=annotation_path_exp,
        combinedWigParams=combinedWigParams,
    ))
    when_data_has_been_collected()

def collect_data_from_gui():
    frame = universal.frame
    with gui_tools.nice_error_log:
        new_parameters = LazyDict(
            combined_wig=None,
            metadata=None,
            annotation=None,
            
            normalization=None,
            output_file=None,
            excluded_conditions=None,
            included_conditions=None,
            n_terminus=None,
            c_terminus=None,
            pseudocount=None,
            winz=None,
            refs=None,
            alpha=None
        )
        
        # 
        # get wig files
        # 
        new_parameters.combined_wig = [ each.cwig.path     for each in universal.session_data.wig_groups ][0]
        new_parameters.metadata     = [ each.metadata.path for each in universal.session_data.wig_groups ][0]
        
        # 
        # get annotation
        # 
        new_parameters.annotation = universal.session_data.annotation
        # TODO: enable this 
        # if not transit_tools.validate_annotation(new_parameters.annotation):
        #     return None
        
        # 
        # setup custom new_parameters
        # 
        
        
        
        # 
        # get output path
        # 
        default_file_name = "test1_output.dat"  # simplified
        default_dir = os.getcwd()
        output_path = frame.SaveFile(default_dir, default_file_name)
        if not output_path:
            return None
        new_parameters.output_file = open(output_path, "w")
        
        
        # 
        # last step
        # 
        parameters.update(new_parameters) # integrate new information
        when_data_has_been_collected() # run

def define_menu_options(himar1_menu, tn5_menu):
    from pytransit.components.menu import method_select_func
    
    def when_menu_clicked(event):
        universal.selected_method = the_method
        # hide all the other panel stuff
        for each_method_name in method_names:
            each_method = methods[each_method_name]
            if each_method.gui.panel:
                each_method.gui.panel.Hide()
        the_method.gui.define_panel(frame)
        return method_select_helper(the_fullname, event)
    
    
    # attach menus
    for method_name, parent_menu in [ ["himar1", himar1_menu_item], ["tn5", tn5_menu_item] ]:
        temp_menu_item = wx.MenuItem(parent_menu, wx.ID_ANY, fullname, wx.EmptyString, wx.ITEM_NORMAL,)
        frame.Bind(wx.EVT_MENU,menu_callback,temp_menu_item,)
        parent_menu.Append(temp_menu_item)


def define_panel():
    from pytransit.components.parameter_panel import panel
    from pytransit.components.panel_helpers import create_normalization_dropdown, create_reference_condition_dropdown, create_include_condition_list
    
    frame = universal.frame
    test1_panel = wx.Panel(
        frame,
        wx.ID_ANY,
        wx.DefaultPosition,
        wx.DefaultSize,
        wx.TAB_TRAVERSAL,
    )
    self.value_getters = LazyDict()
    
    
    # --include-conditions <cond1,...> := Comma-separated list of conditions to use for analysis (Default: all)
    # --exclude-conditions <cond1,...> := Comma-separated list of conditions to exclude (Default: none)
    # --ref <cond> := which condition(s) to use as a reference for calculating LFCs (comma-separated if multiple conditions)
    # -iN <N> :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: -iN 0
    # -iC <N> :=  Ignore TAs within given percentage (e.g. 5) of C terminus. Default: -iC 0
    # -PC <N> := pseudocounts to use for calculating LFC. Default: -PC 5
    # -winz   := winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers)
    
    # 
    # container: parameters
    # 
    if True:
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.value_getters.normalization          = create_normalization_dropdown(test1_panel, main_sizer)
        self.value_getters.reference_condition    = create_reference_condition_dropdown(test1_panel, main_sizer)
        self.value_getters.include_condition_list = create_include_condition_list(test1_panel, main_sizer)
        
        # # 
        # # text input: excludeConditionsInput
        # # 
        # if True:
        #     (
        #         _,
        #         universal.frame.excludeConditionsInput,
        #         sizer,
        #     ) = self.defineTextBox(
        #         panel=test1_panel,
        #         labelText="Exclude\nConditions\n",
        #         widgetText="",
        #         tooltipText="comma seperated list (default=none)",
        #     )
        #     main_sizer.Add(sizer, 1, wx.EXPAND, default_padding)
        
        # # 
        # # text input: Pseudocount
        # # 
        # if True:
        #     # (test1PseudocountLabel, universal.frame.test1PseudocountText, pseudoSizer) = self.defineTextBox(test1_panel, "Pseudocount:", "0.0", "Adds pseudo-counts to the each data-point. Useful to dampen the effects of small counts which may lead to deceptively high log-FC.")
        #     (
        #         test1PseudocountLabel,
        #         universal.frame.test1PseudocountText,
        #         pseudoSizer,
        #     ) = self.defineTextBox(
        #         test1_panel,
        #         "Pseudocount:",
        #         "5",
        #         "Pseudo-counts used in calculating log-fold-change. Useful to dampen the effects of small counts which may lead to deceptively high LFC.",
        #     )
        #     main_sizer.Add(pseudoSizer, 1, wx.EXPAND, default_padding)
        
        # # 
        # # normalization method
        # # 
        # if True:
        #     (
        #         label,
        #         normalization_wxobj,
        #         normalization_choice_sizer,
        #     ) = self.defineChoiceBox(
        #         test1_panel,
        #         "Normalization: ",
        #         [
        #             "TTR",
        #             "nzmean",
        #             "totreads",
        #             "zinfnb",
        #             "quantile",
        #             "betageom",
        #             "nonorm",
        #         ],
        #         "Choice of normalization method. The default choice, 'TTR', normalizes datasets to have the same expected count (while not being sensative to outliers). Read documentation for a description other methods. ",
        #     )
        #     main_sizer.Add(normalization_choice_sizer, 1, wx.EXPAND, default_padding)
        #     self.get_normalization_choice = lambda : normalization_wxobj.GetString(normalization_wxobj.GetCurrentSelection())

        # main_sizer.Add(main_sizer, 1, wx.EXPAND, default_padding)
    
    # 
    # checkbox: Correct for Genome Positional Bias
    # 
    # if True:
        # LOESS Check
        # (universal.frame.test1LoessCheck, loessCheckSizer) = self.defineCheckBox(
        #     test1_panel,
        #     labelText="Correct for Genome Positional Bias",
        #     widgetCheck=False,
        #     widgetSize=(-1, -1),
        #     tooltipText="Check to correct read-counts for possible regional biase using LOESS. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.",
        # )
        # main_sizer.Add(loessCheckSizer, 0, wx.EXPAND, default_padding)
    
    # 
    # button: LOESS
    # 
    # if True:
    #     universal.frame.test1LoessPrev = wx.Button(
    #         test1_panel,
    #         wx.ID_ANY,
    #         "Preview LOESS fit",
    #         wx.DefaultPosition,
    #         wx.DefaultSize,
    #         0,
    #     )
    #     universal.frame.test1LoessPrev.Bind(wx.EVT_BUTTON, universal.frame.LoessPrevFunc)
    
    # # 
    # # checkbox: winsorize
    # # 
    # if True:
    #     (universal.frame.test1AdaptiveCheckBox, adaptiveSizer) = self.defineCheckBox(
    #         test1_panel,
    #         labelText="winsorize",
    #         widgetCheck=False,
    #         widgetSize=(-1, -1),
    #         tooltipText="winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers)",
    #     )
    #     main_sizer.Add(adaptiveSizer, 1, wx.ALL | wx.EXPAND, default_padding)
        
        # 
        # button: RUN
        # 
        if True:
            runButton = wx.Button(
                test1_panel,
                wx.ID_ANY,
                "Run",
                wx.DefaultPosition,
                wx.DefaultSize,
                0,
            )
            @gui_tools.bind_to(runButton, wx.EVT_BUTTON)
            def when_clicked(*args):
                with gui_tools.nice_error_log:
                    # pull in the GUI data
                    from_gui()
                    # run the process in the background
                    thread = threading.Thread(target=Run())
                    thread.setDaemon(True)
                    thread.start()
            main_sizer.Add(runButton, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, default_padding)
        
        panel.method_sizer.Add(test1_panel, 0, wx.EXPAND, default_padding)
        
        test1_panel.SetSizer(main_sizer)
        test1_panel.Layout()
        main_sizer.Fit(test1_panel)

    return test1_panel


# 
# 
# Helpers
# 
# 

def preprocess_data(position, data):
    (K, N) = data.shape

    if parameters.normalization != "nonorm":
        gui_tools.set_status("Normalizing using: %s" % parameters.normalization)
        (data, factors) = norm_tools.normalize_data(
            data,
            parameters.normalization,
            parameters.ctrldata + parameters.expdata,
            parameters.annotation_path,
        )

    if parameters.LOESS:
        gui_tools.set_status("Performing LOESS Correction")
        for j in range(K):
            data[j] = stat_tools.loess_correction(position, data[j])

    return data

def wigs_to_conditions(conditionsByFile, filenamesInCombWig):
    """
        Returns list of conditions corresponding to given wigfiles.
        ({FileName: Condition}, [FileName]) -> [Condition]
        Condition :: [String]
    """
    return [conditionsByFile.get(f, None) for f in filenamesInCombWig]

def filter_wigs_by_conditions(data, conditions, included_conditions):
    """
        Filters conditions from wig to ctrl, exp conditions only
        ([[Wigdata]], [ConditionCtrl, ConditionExp]) -> Tuple([[Wigdata]], [Condition])
    """
    d_filtered, cond_filtered = [], []
    if len(included_conditions) != 2:
        raise Exception(f'''Only 2 conditions expected: {included_conditions}''')
    for i, c in enumerate(conditions):
        if c.lower() in included_conditions:
            d_filtered.append(data[i])
            cond_filtered.append(conditions[i])

    return (numpy.array(d_filtered), numpy.array(cond_filtered))

def winsorize_test1(self, counts):
    # input is insertion counts for gene as pre-flattened numpy array
    counts = counts.tolist()
    if len(counts) < 3:
        return counts
    s = sorted(counts, reverse=True)
    if s[1] == 0:
        return counts  # don't do anything if there is only 1 non-zero value
    c2 = [s[1] if x == s[0] else x for x in counts]
    return numpy.array(c2)

    # unique_counts = numpy.unique(counts)
    # if (len(unique_counts) < 2): return counts
    # else:
    #  n, n_minus_1 = unique_counts[heapq.nlargest(2, range(len(unique_counts)), unique_counts.take)]
    #  result = [[ n_minus_1 if count == n else count for count in wig] for wig in counts]
    #  return numpy.array(result)

def run_test1(
    self, G_ctrl, G_exp=None, doLibraryAnova=False, histPath=""
):
    data = []
    N = len(G_ctrl)
    count = 0
    self.progress_range(N)

    for gene in G_ctrl:
        if gene.orf not in G_exp:
            if self.diffStrains:
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

        if not self.diffStrains and gene.n != gene_exp.n:
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
            if not self.includeZeros:
                ii_ctrl = numpy.sum(gene.reads, 0) > 0
                ii_exp = numpy.sum(gene_exp.reads, 0) > 0
            else:
                ii_ctrl = numpy.ones(gene.n) == 1
                ii_exp = numpy.ones(gene_exp.n) == 1

            # data1 = gene.reads[:,ii_ctrl].flatten() + self.pseudocount # we used to have an option to add pseudocounts to each observation, like this
            data1 = gene.reads[:, ii_ctrl].flatten()
            data2 = gene_exp.reads[:, ii_exp].flatten()
            if self.winz:
                data1 = self.winsorize_test1(data1)
                data2 = self.winsorize_test1(data2)

            if doLibraryAnova:
                (
                    test_obs,
                    mean1,
                    mean2,
                    log2FC,
                    pval_ltail,
                    pval_utail,
                    pval_2tail,
                    testlist,
                ) = stat_tools.test1(
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
                ) = stat_tools.test1(
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

        if self.doHistogram:
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
        text = "Running Anova Method... %5.1f%%" % (100.0 * count / N)
        self.progress_update(text, count)

    gui_tools.set_status("")  # Printing empty line to flush stdout
    gui_tools.set_status("Performing Benjamini-Hochberg Correction")
    data.sort()
    qval = stat_tools.BH_fdr_correction([row[-1] for row in data])

    return (data, qval)


self.fullname = f"[{self.short_name}]  -  {self.short_desc}"