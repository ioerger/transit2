import sys

from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub

import os
import time
import ntpath
import math
import random
import numpy
import scipy.stats
import datetime
import heapq

from pytransit.analysis import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools


############# GUI ELEMENTS ##################

short_name = "test1"
long_name = "Test1 (Permutation test)"
short_desc = "Test1 test of conditional essentiality between two conditions"
long_desc = """Method for determining conditional essentiality based on test1 (i.e. permutation test). Identifies significant changes in mean read-counts for each gene after normalization."""

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


class Analysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            transposons,
            Method,
            GUI,
            [File],
        )


############# FILE ##################


class File(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#Test1", columns)

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
            transit_tools.ShowError(
                MSG="Error Displaying File. Histogram image not found. Make sure results were obtained with the histogram option turned on."
            )
            print("Error Displaying File. Histogram image does not exist.")


############# GUI ##################


class GUI(base.AnalysisGUI):
    def definePanel(self, wxobj):
        self.wxobj = wxobj
        test1Panel = wx.Panel(
            self.wxobj.optionsWindow,
            wx.ID_ANY,
            wx.DefaultPosition,
            wx.DefaultSize,
            wx.TAB_TRAVERSAL,
        )
        
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
            test1Sizer = wx.BoxSizer(wx.VERTICAL)
            
            # 
            # label
            # 
            if True:
                test1Label = wx.StaticText(
                    test1Panel,
                    wx.ID_ANY,
                    u"Local Options",
                    wx.DefaultPosition,
                    (160, -1),
                    0,
                )
                test1Label.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
                test1Sizer.Add(test1Label, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
            
            # 
            # main sizer
            # 
            if True:
                mainSizer1 = wx.BoxSizer(wx.VERTICAL)
                
                # 
                # text input: includeConditionsInput
                # 
                if True:
                    (
                        _,
                        self.wxobj.includeConditionsInput,
                        sizer,
                    ) = self.defineTextBox(
                        panel=test1Panel,
                        labelText=u"Include\nConditions\n",
                        widgetText=u"",
                        tooltipText="comma seperated list (default=all)",
                    )
                    mainSizer1.Add(sizer, 1, wx.EXPAND, 5)

                # 
                # text input: excludeConditionsInput
                # 
                if True:
                    (
                        _,
                        self.wxobj.excludeConditionsInput,
                        sizer,
                    ) = self.defineTextBox(
                        panel=test1Panel,
                        labelText=u"Exclude\nConditions\n",
                        widgetText=u"",
                        tooltipText="comma seperated list (default=none)",
                    )
                    mainSizer1.Add(sizer, 1, wx.EXPAND, 5)
                
                # 
                # text input: refInput
                # 
                if True:
                    # TODO: could be dropdown
                    (
                        _,
                        self.wxobj.refInput,
                        sizer,
                    ) = self.defineTextBox(
                        panel=test1Panel,
                        labelText=u"Ref",
                        widgetText=u"",
                        tooltipText="which condition(s) to use as a reference for calculating LFCs (comma-separated if multiple conditions)",
                    )
                    mainSizer1.Add(sizer, 1, wx.EXPAND, 5)


                # 
                # text input: Pseudocount
                # 
                if True:
                    # (test1PseudocountLabel, self.wxobj.test1PseudocountText, pseudoSizer) = self.defineTextBox(test1Panel, u"Pseudocount:", u"0.0", "Adds pseudo-counts to the each data-point. Useful to dampen the effects of small counts which may lead to deceptively high log-FC.")
                    (
                        test1PseudocountLabel,
                        self.wxobj.test1PseudocountText,
                        pseudoSizer,
                    ) = self.defineTextBox(
                        test1Panel,
                        u"Pseudocount:",
                        u"0.0",
                        "Pseudo-counts used in calculating log-fold-chnage. Useful to dampen the effects of small counts which may lead to deceptively high LFC.",
                    )
                    mainSizer1.Add(pseudoSizer, 1, wx.EXPAND, 5)
                
                # 
                # text input: 
                # 
                if True:
                    # -iN <N> :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: -iN 0
                    (
                        _,
                        self.wxobj.nTerminusInput,
                        sizer,
                    ) = self.defineTextBox(
                        test1Panel,
                        u"N Terminus\nTruncate",
                        u"0.0",
                        "Ignore TAs within given percentage (e.g. 5) of N terminus. Default: 0",
                    )
                    mainSizer1.Add(sizer, 1, wx.EXPAND, 5)
                
                # 
                # text input: 
                # 
                if True:
                    (
                        _,
                        self.wxobj.cTerminusInput,
                        sizer,
                    ) = self.defineTextBox(
                        test1Panel,
                        u"C Terminus\nTruncate",
                        u"0.0",
                        "Ignore TAs within given percentage (e.g. 5) of C terminus. Default: 0",
                    )
                    mainSizer1.Add(sizer, 1, wx.EXPAND, 5)
                
                # 
                # dropdown
                # 
                if True:
                    test1NormChoiceChoices = [
                        u"TTR",
                        u"nzmean",
                        u"totreads",
                        u"zinfnb",
                        u"quantile",
                        u"betageom",
                        u"nonorm",
                    ]
                    (
                        test1NormLabel,
                        self.wxobj.test1NormChoice,
                        normSizer,
                    ) = self.defineChoiceBox(
                        test1Panel,
                        u"Normalization: ",
                        test1NormChoiceChoices,
                        "Choice of normalization method. The default choice, 'TTR', normalizes datasets to have the same expected count (while not being sensative to outliers). Read documentation for a description other methods. ",
                    )
                    mainSizer1.Add(normSizer, 1, wx.EXPAND, 5)

                test1Sizer.Add(mainSizer1, 1, wx.EXPAND, 5)
            
            # 
            # checkbox: Correct for Genome Positional Bias
            # 
            # if True:
                # LOESS Check
                # (self.wxobj.test1LoessCheck, loessCheckSizer) = self.defineCheckBox(
                #     test1Panel,
                #     labelText="Correct for Genome Positional Bias",
                #     widgetCheck=False,
                #     widgetSize=(-1, -1),
                #     tooltipText="Check to correct read-counts for possible regional biase using LOESS. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.",
                # )
                # test1Sizer.Add(loessCheckSizer, 0, wx.EXPAND, 5)
            
            # 
            # button: LOESS
            # 
            # if True:
            #     self.wxobj.test1LoessPrev = wx.Button(
            #         test1Panel,
            #         wx.ID_ANY,
            #         u"Preview LOESS fit",
            #         wx.DefaultPosition,
            #         wx.DefaultSize,
            #         0,
            #     )
            #     self.wxobj.test1LoessPrev.Bind(wx.EVT_BUTTON, self.wxobj.LoessPrevFunc)
            
            # 
            # checkbox: winsorize
            # 
            if True:
                (self.wxobj.test1AdaptiveCheckBox, adaptiveSizer) = self.defineCheckBox(
                    test1Panel,
                    labelText="winsorize",
                    widgetCheck=False,
                    widgetSize=(-1, -1),
                    tooltipText="winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers)",
                )
                test1Sizer.Add(adaptiveSizer, 0, wx.EXPAND, 5)
            
            # 
            # button: RUN
            # 
            if True:
                runButton = wx.Button(
                    test1Panel,
                    wx.ID_ANY,
                    u"Run",
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                # runButton.SetBackgroundColour((0, 230, 200, 255)) # rgba(0, 230, 200, 255)
                # runButton.SetForegroundColour((70, 70, 70, 255)) # rgba(70, 70, 70, 255)
                runButton.Bind(wx.EVT_BUTTON, self.wxobj.RunMethod)
                test1Sizer.Add(runButton, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

            test1Panel.SetSizer(test1Sizer)
            test1Panel.Layout()
            test1Sizer.Fit(test1Panel)

        self.panel = test1Panel

    def GlobalEnable(self):
        self.wxobj.ctrlLibText.Enable()
        self.wxobj.expLibText.Enable()

    def GlobalDisable(self):
        self.wxobj.ctrlLibText.Disable()
        self.wxobj.expLibText.Disable()


########## CLASS #######################


class Method(base.DualConditionMethod):
    """
    test1

    """

    def __init__(
        self,
        ctrldata,
        expdata,
        annotation_path,
        output_file,
        normalization="TTR",
        samples=10000,
        adaptive=False,
        doHistogram=False,
        includeZeros=False,
        pseudocount=1,
        replicates="Sum",
        LOESS=False,
        ignoreCodon=True,
        NTerminus=0.0,
        CTerminus=0.0,
        ctrl_lib_str="",
        exp_lib_str="",
        winz=False,
        wxobj=None,
        Z=False,
        diffStrains=False,
        annotation_path_exp="",
        combinedWigParams=None,
    ):

        base.DualConditionMethod.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            ctrldata,
            expdata,
            annotation_path,
            output_file,
            normalization=normalization,
            replicates=replicates,
            LOESS=LOESS,
            NTerminus=NTerminus,
            CTerminus=CTerminus,
            wxobj=wxobj,
        )

        self.Z = Z
        self.samples = samples
        self.adaptive = adaptive
        self.doHistogram = doHistogram
        self.includeZeros = includeZeros
        self.pseudocount = pseudocount
        self.ctrl_lib_str = ctrl_lib_str
        self.exp_lib_str = exp_lib_str
        self.diffStrains = diffStrains
        self.annotation_path_exp = (
            annotation_path_exp if diffStrains else annotation_path
        )
        self.combinedWigParams = combinedWigParams
        self.winz = winz

    @classmethod
    def fromGUI(self, wxobj):
        """ """
        # Get Annotation file
        annot_paths = wxobj.annotation.split(",")
        annotationPath = annot_paths[0]
        diffStrains = False
        annotationPathExp = ""
        if len(annot_paths) == 2:
            annotationPathExp = annot_paths[1]
            diffStrains = True

        if not transit_tools.validate_annotation(annotationPath):
            return None

        if annotationPathExp and not transit_tools.validate_annotation(
            annotationPathExp
        ):
            return None

        # Get selected files
        ctrldata = wxobj.ctrlSelected()
        expdata = wxobj.expSelected()
        if not transit_tools.validate_both_datasets(ctrldata, expdata):
            return None

        # Validate transposon types
        if not transit_tools.validate_transposons_used(ctrldata + expdata, transposons):
            return None

        # Read the parameters from the wxPython widgets
        ignoreCodon = True
        samples = int(wxobj.normalizationMethodInput.GetValue())
        normalization = wxobj.test1NormChoice.GetString(
            wxobj.test1NormChoice.GetCurrentSelection()
        )
        replicates = "Sum"
        adaptive = wxobj.test1AdaptiveCheckBox.GetValue()
        doHistogram = wxobj.test1HistogramCheckBox.GetValue()

        includeZeros = wxobj.test1ZeroCheckBox.GetValue()
        pseudocount = float(wxobj.test1PseudocountText.GetValue())
        LOESS = wxobj.test1LoessCheck.GetValue()

        # Global Parameters
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        ctrl_lib_str = wxobj.ctrlLibText.GetValue()
        exp_lib_str = wxobj.expLibText.GetValue()

        # Get output path
        # defaultFileName = "test1_output_s%d_pc%1.2f" % (samples, pseudocount)
        # if adaptive: defaultFileName+= "_adaptive"
        # if includeZeros: defaultFileName+= "_iz"
        # defaultFileName+=".dat"
        defaultFileName = "test1_output.dat"  # simplified

        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path:
            return None
        output_file = open(output_path, "w")

        return self(
            ctrldata,
            expdata,
            annotationPath,
            output_file,
            normalization,
            samples,
            adaptive,
            doHistogram,
            includeZeros,
            pseudocount,
            replicates,
            LOESS,
            ignoreCodon,
            NTerminus,
            CTerminus,
            ctrl_lib_str,
            exp_lib_str,
            wxobj,
            Z=False,
            diffStrains=diffStrains,
            annotation_path_exp=annotationPathExp,
        )

    @classmethod
    def fromargs(self, rawargs):

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
        annotationPath = annot_paths[0]
        diffStrains = False
        annotationPathExp = ""
        if len(annot_paths) == 2:
            annotationPathExp = annot_paths[1]
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
        samples = int(kwargs.get("s", 10000))
        adaptive = kwargs.get("a", False)
        doHistogram = kwargs.get("h", False)
        replicates = kwargs.get("r", "Sum")
        excludeZeros = kwargs.get("ez", False)
        includeZeros = not excludeZeros
        pseudocount = float(
            kwargs.get("PC", 1.0)
        )  # use -PC (new semantics: for LFCs) instead of -pc (old semantics: fake counts)

        Z = True if "Z" in kwargs else False

        LOESS = kwargs.get("l", False)
        ignoreCodon = True

        NTerminus = float(kwargs.get("iN", 0.00))  # integer interpreted as percentage
        CTerminus = float(kwargs.get("iC", 0.00))
        ctrl_lib_str = kwargs.get("-ctrl_lib", "")
        exp_lib_str = kwargs.get("-exp_lib", "")

        return self(
            ctrldata,
            expdata,
            annotationPath,
            output_file,
            normalization,
            samples,
            adaptive,
            doHistogram,
            includeZeros,
            pseudocount,
            replicates,
            LOESS,
            ignoreCodon,
            NTerminus,
            CTerminus,
            ctrl_lib_str,
            exp_lib_str,
            winz=winz,
            Z=Z,
            diffStrains=diffStrains,
            annotation_path_exp=annotationPathExp,
            combinedWigParams=combinedWigParams,
        )

    def preprocess_data(self, position, data):
        (K, N) = data.shape

        if self.normalization != "nonorm":
            self.transit_message("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(
                data,
                self.normalization,
                self.ctrldata + self.expdata,
                self.annotation_path,
            )

        if self.LOESS:
            self.transit_message("Performing LOESS Correction")
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

    def Run(self):

        # if not self.wxobj:
        #    # Force matplotlib to use good backend for png.
        #    import matplotlib.pyplot as plt
        # elif "matplotlib.pyplot" not in sys.modules:
        try:
            import matplotlib.pyplot as plt
        except:
            print("Error: cannot do histograms")
            self.doHistogram = False

        self.transit_message("Starting test1 Method")
        start_time = time.time()
        if self.winz:
            self.transit_message("Winsorizing insertion counts")

        histPath = ""
        if self.doHistogram:
            histPath = os.path.join(
                os.path.dirname(self.output.name),
                transit_tools.fetch_name(self.output.name) + "_histograms",
            )
            if not os.path.isdir(histPath):
                os.makedirs(histPath)

        # Get orf data
        self.transit_message("Getting Data")
        if self.diffStrains:
            self.transit_message("Multiple annotation files found")
            self.transit_message(
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
                self.ctrldata, wxobj=self.wxobj
            )
            (data_exp, position_exp) = transit_tools.get_validated_data(
                self.expdata, wxobj=self.wxobj
            )
        (K_ctrl, N_ctrl) = data_ctrl.shape
        (K_exp, N_exp) = data_exp.shape

        if not self.diffStrains and (N_ctrl != N_exp):
            self.transit_error(
                "Error: Ctrl and Exp wig files don't have the same number of sites."
            )
            self.transit_error("Make sure all .wig files come from the same strain.")
            return
        # (data, position) = transit_tools.get_validated_data(self.ctrldata+self.expdata, wxobj=self.wxobj)

        self.transit_message("Preprocessing Ctrl data...")
        data_ctrl = self.preprocess_data(position_ctrl, data_ctrl)

        self.transit_message("Preprocessing Exp data...")
        data_exp = self.preprocess_data(position_exp, data_exp)

        G_ctrl = tnseq_tools.Genes(
            self.ctrldata,
            self.annotation_path,
            ignoreCodon=self.ignoreCodon,
            nterm=self.NTerminus,
            cterm=self.CTerminus,
            data=data_ctrl,
            position=position_ctrl,
        )
        G_exp = tnseq_tools.Genes(
            self.expdata,
            self.annotation_path_exp,
            ignoreCodon=self.ignoreCodon,
            nterm=self.NTerminus,
            cterm=self.CTerminus,
            data=data_exp,
            position=position_exp,
        )

        doLibraryTest1 = False
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
                    doLibraryTest1 = True
                else:
                    transit_tools.transit_error(
                        "Error: Library Strings (Ctrl = %s, Exp = %s) do not use the same letters. Make sure every letter / library is represented in both Control and Experimental Conditions. Proceeding with test1 assuming all datasets belong to the same library."
                        % (self.ctrl_lib_str, self.exp_lib_str)
                    )
                    self.ctrl_lib_str = ""
                    self.exp_lib_str = ""

        (data, qval) = self.run_test1(G_ctrl, G_exp, doLibraryTest1, histPath)
        self.write_output(data, qval, start_time)

        self.finish()
        self.transit_message("Finished test1 Method")

    def write_output(self, data, qval, start_time):

        self.output.write("#Test1\n")
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
                "#GUI with: norm=%s, samples=%s, pseudocounts=%1.2f, adaptive=%s, histogram=%s, includeZeros=%s, output=%s\n"
                % (
                    self.normalization,
                    self.samples,
                    self.pseudocount,
                    self.adaptive,
                    self.doHistogram,
                    self.includeZeros,
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
                self.doHistogram,
                self.adaptive,
                not self.includeZeros,
                self.pseudocount,
                self.LOESS,
                self.NTerminus,
                self.CTerminus,
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
                self.annotation_path_exp.encode("utf-8") if self.diffStrains else "",
            )
        )
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        # Z = True # include Z-score column in test1 output?
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
                    p = 1e-5  # or 1 level deeper the num of iterations of test1, which is 1e-4=1/10000, by default
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

        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="Test1")

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
        self, G_ctrl, G_exp=None, doLibraryTest1=False, histPath=""
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

                if doLibraryTest1:
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
                        PC=self.pseudocount,
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
                        PC=self.pseudocount,
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
            text = "Running Test1 Method... %5.1f%%" % (100.0 * count / N)
            self.progress_update(text, count)

        #
        self.transit_message("")  # Printing empty line to flush stdout
        self.transit_message("Performing Benjamini-Hochberg Correction")
        data.sort()
        qval = stat_tools.BH_fdr_correction([row[-1] for row in data])

        return (data, qval)

    @classmethod
    def usage_string(self):
        return """
        python3 %s test1 <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        ---
        OR
        ---
        python3 %s test1 -c <combined wig file> <samples_metadata file> <ctrl condition name> <exp condition name> <annotation .prot_table> <output file> [Optional Arguments]
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
        """ % (
            sys.argv[0],
            sys.argv[0],
        )


if __name__ == "__main__":

    (args, kwargs) = transit_tools.clean_args(sys.argv)

    # TODO: Figure out issue with inputs (transit requires initial method name, running as script does not !!!!)

    G = Method.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")
    G.print_members()

    print("")
    print("Running:")

    G.Run()