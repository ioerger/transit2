import scipy
import numpy
import heapq
import math
import statsmodels.stats.multitest

import time
import sys
import collections

from pytransit.analysis import base
from pytransit.transit_tools import write_dat, EOL
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools


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
import pytransit.components.results_area as results_area
from pytransit.core_data import universal
from pytransit.components.parameter_panel import panel
from pytransit.components.panel_helpers import make_panel, create_run_button, create_normalization_input, create_reference_condition_input, create_include_condition_list_input, create_exclude_condition_list_input, create_n_terminus_input, create_c_terminus_input, create_pseudocount_input, create_winsorize_input, create_alpha_input

############# GUI ELEMENTS ##################

main_object = LazyDict(
    short_name = "anova -- test",
    long_name = "AnovaGUI",
    short_desc = "AnovaGUI",
    long_desc = """Anova GUI""",

    transposons = ["himar1", "tn5"],
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
    ],
    
    analysis=None,
    gui=None,
    method=None,
    
    inputs = LazyDict(
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
        alpha=None,
        combinedWigParams=None,
    ),
)

class Analysis(base.TransitAnalysis):
    def __init__(self):
        main_object.analysis = self
        base.TransitAnalysis.__init__(
            self,
            
            main_object.short_name,
            main_object.long_name,
            main_object.short_desc,
            main_object.long_desc,
            main_object.transposons,
            
            Method,
            GUI,
            [File],
        )

class File(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#Anova", columns)

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

class GUI:
    def __init__(self, *args, **kwargs):
        main_object.gui = self
        self.wxobj = None
        self.panel = None
    
    def define_panel(self, _):
        
        self.main_frame = frame = universal.frame
        test1_panel = make_panel()
        
        # 
        # parameter inputs
        # 
        # --include-conditions <cond1,...> := Comma-separated list of conditions to use for analysis (Default: all)
        # --exclude-conditions <cond1,...> := Comma-separated list of conditions to exclude (Default: none)
        # --ref <cond> := which condition(s) to use as a reference for calculating LFCs (comma-separated if multiple conditions)
        # -iN <N> :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: -iN 0
        # -iC <N> :=  Ignore TAs within given percentage (e.g. 5) of C terminus. Default: -iC 0
        # -PC <N> := pseudocounts to use for calculating LFC. Default: -PC 5
        # -winz   := winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers)
        self.value_getters = LazyDict()
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        if True:
            self.value_getters.included_conditions    = create_include_condition_list_input(test1_panel, main_sizer)
            self.value_getters.excluded_conditions    = create_exclude_condition_list_input(test1_panel, main_sizer)
            self.value_getters.reference_condition    = create_reference_condition_input(test1_panel, main_sizer)
            self.value_getters.n_terminus             = create_n_terminus_input(test1_panel, main_sizer)
            self.value_getters.c_terminus             = create_c_terminus_input(test1_panel, main_sizer)
            self.value_getters.normalization          = create_normalization_input(test1_panel, main_sizer)
            self.value_getters.pseudocount            = create_pseudocount_input(test1_panel, main_sizer)
            self.value_getters.alpha                  = create_alpha_input(test1_panel, main_sizer)
            self.value_getters.winz                   = create_winsorize_input(test1_panel, main_sizer)
            self.value_getters.refs                   = lambda *args: [] if self.value_getters.reference_condition() == "[None]" else [ self.value_getters.reference_condition() ]
            
            create_run_button(test1_panel, main_sizer)
            
        panel.method_sizer.Add(test1_panel, 0, wx.EXPAND, gui_tools.default_padding)
        test1_panel.SetSizer(main_sizer)
        test1_panel.Layout()
        main_sizer.Fit(test1_panel)

        self.panel = test1_panel

class Method(base.MultiConditionMethod):
    def __init__(
        self,
        combined_wig,
        metadata,
        annotation,
        normalization,
        output_file,
        excluded_conditions=[],
        included_conditions=[],
        n_terminus=0.0,
        c_terminus=0.0,
        pseudocount=1,
        winz=False,
        refs=[],
        alpha=1000,
        *args,
        **kwargs,
    ):
        main_object.method = self
        
        self.pseudocount         = None
        self.alpha               = None
        self.refs                = None
        self.winz                = None
        self.normalization       = None
        self.ctrldata            = None
        self.expdata             = None
        self.annotation_path     = None
        self.LOESS               = None
        self.main_frame          = None
        self.doHistogram         = None
        self.output              = None
        self.diffStrains         = None
        self.annotation_path_exp = None
        self.combinedWigParams   = None
        self.ignoreCodon         = None
        self.n_terminus          = None
        self.c_terminus          = None
        self.ctrl_lib_str        = None
        self.exp_lib_str         = None
        self.samples             = None
        self.adaptive            = None
        self.includeZeros        = None
        self.Z                   = None
        
        base.MultiConditionMethod.__init__(
            self,
            main_object.short_name,
            main_object.long_name,
            main_object.short_desc,
            main_object.long_desc,
            combined_wig,
            metadata,
            annotation,
            output_file,
            normalization=normalization,
            excluded_conditions=excluded_conditions,
            included_conditions=included_conditions,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
        )

        self.pseudocount = pseudocount
        self.alpha = alpha
        self.refs = refs
        self.winz = winz
    
    def __repr__(self):
        as_dict = LazyDict(
            pseudocount=self.pseudocount,
            alpha=self.alpha,
            refs=self.refs,
            winz=self.winz,
            normalization=self.normalization,
            ctrldata=self.ctrldata,
            expdata=self.expdata,
            annotation_path=self.annotation_path,
            LOESS=self.LOESS,
            main_frame=self.main_frame,
            doHistogram=self.doHistogram,
            output=self.output,
            diffStrains=self.diffStrains,
            annotation_path_exp=self.annotation_path_exp,
            combinedWigParams=self.combinedWigParams,
            ignoreCodon=self.ignoreCodon,
            n_terminus=self.n_terminus,
            c_terminus=self.c_terminus,
            ctrl_lib_str=self.ctrl_lib_str,
            exp_lib_str=self.exp_lib_str,
            samples=self.samples,
            adaptive=self.adaptive,
            includeZeros=self.includeZeros,
            Z=self.Z,
        )
        return f"{as_dict}"

    @classmethod
    def usage_string(cls):
        command_name = sys.argv[0]
        return f"""
        python3 {command_name} test1 <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        ---
        OR
        ---
        python3 {command_name} test1 -c <combined wig file> <samples_metadata file> <ctrl condition name> <exp condition name> <annotation .prot_table> <output file> [Optional Arguments]
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
        """

    
    @classmethod
    def from_gui(self, frame):
        with gui_tools.nice_error_log:
            # 
            # get wig files
            # 
            wig_group = universal.session_data.wig_groups[0]
            main_object.inputs.combined_wig = wig_group.cwig.path
            main_object.inputs.metadata     = wig_group.metadata.path
            
            # 
            # get annotation
            # 
            main_object.inputs.annotation = universal.session_data.annotation
            # FIXME: enable this once I get a valid annotation file example
            # if not transit_tools.validate_annotation(main_object.inputs.annotation):
            #     return None
            
            # 
            # setup custom inputs
            # 
            for each_key, each_getter in main_object.gui.value_getters.items():
                try:
                    main_object.inputs[each_key] = each_getter()
                except Exception as error:
                    raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
            
            # 
            # save result files
            # 
            main_object.inputs.output_file = gui_tools.ask_for_output_file_path(
                default_file_name="test1_output.dat",
                output_extensions=u'Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
            )
            if not main_object.inputs.output_file:
                return None

            return self(*main_object.inputs.values())

    @classmethod
    def fromargs(self, rawargs):
        (args, kwargs) = transit_tools.clean_args(rawargs)

        if kwargs.get("-help", False) or kwargs.get("h", False):
            print(self.usage_string)
            sys.exit(0)

        combined_wig = args[0]
        annotation   = args[2]
        metadata     = args[1]
        output_file  = args[3]
        normalization = kwargs.get("n", "TTR")
        n_terminus    = float(kwargs.get("iN", 0.0))
        c_terminus    = float(kwargs.get("iC", 0.0))
        winz          = "winz" in kwargs
        pseudocount   = int(kwargs.get("PC", 5))
        alpha         = float(kwargs.get("alpha", 1000))
        refs          = kwargs.get("-ref", [])  # list of condition names to use a reference for calculating LFCs
        if refs != []: refs = refs.split(",")
        excluded_conditions = list( filter(None, kwargs.get("-exclude-conditions", "").split(",")) )
        included_conditions = list( filter(None, kwargs.get("-include-conditions", "").split(",")) )

        # check for unrecognized flags
        flags = "-n --exclude-conditions --include-conditions -iN -iC -PC --ref -winz -alpha".split()
        for arg in rawargs:
            if arg[0] == "-" and arg not in flags:
                self.transit_error("flag unrecognized: %s" % arg)
                print(self.usage_string)
                sys.exit(0)

        return self(
            combined_wig,
            metadata,
            annotation,
            normalization,
            output_file,
            excluded_conditions,
            included_conditions,
            n_terminus,
            c_terminus,
            pseudocount,
            winz,
            refs,
            alpha
        )

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

        self.output.write("#Anova\n")
        if self.main_frame:
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

        transit_tools.log("Adding File: %s" % (self.output.name))
        self.add_file(filetype="Anova")

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

    def means_by_condition_for_gene(self, sites, conditions, data):
        """
            Returns a dictionary of {Condition: Mean} for each condition.
            ([Site], [Condition]) -> {Condition: Number}
            Site :: Number
            Condition :: String
        """
        nTASites = len(sites)
        wigsByConditions = collections.defaultdict(lambda: [])
        for i, c in enumerate(conditions):
            wigsByConditions[c].append(i)

        return {
            c: numpy.mean(self.winsorize(data[wigIndex][:, sites]))
            if nTASites > 0
            else 0
            for (c, wigIndex) in wigsByConditions.items()
        }

    def means_by_rv(self, data, RvSiteindexesMap, genes, conditions):
        """
            Returns Dictionary of mean values by condition
            ([[Wigdata]], {Rv: SiteIndex}, [Gene], [Condition]) -> {Rv: {Condition: Number}}
            Wigdata :: [Number]
            SiteIndex :: Number
            Gene :: {start, end, rv, gene, strand}
            Condition :: String
        """
        MeansByRv = {}
        for gene in genes:
            Rv = gene["rv"]
            MeansByRv[Rv] = self.means_by_condition_for_gene(
                RvSiteindexesMap[Rv], conditions, data
            )
        return MeansByRv

    def group_by_condition(self, wigList, conditions):
        """
            Returns array of datasets, where each dataset corresponds to one condition.
            ([[Wigdata]], [Condition]) -> [[DataForCondition]]
            Wigdata :: [Number]
            Condition :: String
            DataForCondition :: [Number]
        """
        countsByCondition = collections.defaultdict(lambda: [])
        countSum = 0
        for i, c in enumerate(conditions):
            countSum += numpy.sum(wigList[i])
            countsByCondition[c].append(wigList[i])

        return (
            countSum,
            [numpy.array(v).flatten() for v in countsByCondition.values()],
        )

    # since this is in both ZINB and ANOVA, should move to stat_tools.py?

    def winsorize(self, counts):
        # input is insertion counts for gene: list of lists: n_replicates (rows) X n_TA sites (cols) in gene
        unique_counts = numpy.unique(numpy.concatenate(counts))
        if len(unique_counts) < 2:
            return counts
        else:
            n, n_minus_1 = unique_counts[
                heapq.nlargest(2, range(len(unique_counts)), unique_counts.take)
            ]
            result = [
                [n_minus_1 if count == n else count for count in wig] for wig in counts
            ]
            return numpy.array(result)

    def run_anova(self, data, genes, MeansByRv, RvSiteindexesMap, conditions):
        """
            Runs Anova (grouping data by condition) and returns p and q values
            ([[Wigdata]], [Gene], {Rv: {Condition: Mean}}, {Rv: [SiteIndex]}, [Condition]) -> Tuple([Number], [Number])
            Wigdata :: [Number]
            Gene :: {start, end, rv, gene, strand}
            Mean :: Number
            SiteIndex: Integer
            Condition :: String
        """
        count = 0
        self.progress_range(len(genes))

        MSR, MSE, Fstats, pvals, Rvs, status = [],[],[],[],[],[]
        for gene in genes:
            count += 1
            Rv = gene["rv"]
            if len(RvSiteindexesMap[Rv]) <= 1:
                status.append("TA sites <= 1")
                msr,mse,Fstat,pval = 0,0,-1,1
            else:
                countSum, countsVec = self.group_by_condition(
                    list(map(lambda wigData: wigData[RvSiteindexesMap[Rv]], data)),
                    conditions,
                )
                if self.winz:
                    countsVec = self.winsorize(countsVec)

                if countSum == 0:
                    msr,mse,Fstat,pval = 0,0,-1,1
                    status.append("No counts in all conditions")
                else:
                    Fstat,pval = scipy.stats.f_oneway(*countsVec)
                    status.append("-")
                    # countsVec is a list of numpy arrays, or could be a list of lists
                    # pooled counts for each condition, over TAs in gene and replicates
                    if isinstance(countsVec[0],numpy.ndarray): 
                      countsVecAsArrays = countsVec
                      countsVecAsLists = [grp.tolist() for grp in countsVec]
                    else:
                      countsVecAsArrays = [numpy.array(grp) for grp in countsVec]
                      countsVecAsLists = countsVec
                    allcounts = [item for sublist in countsVecAsLists for item in sublist]
                    grandmean = numpy.mean(allcounts)
                    groupmeans = [numpy.mean(grp) for grp in countsVecAsArrays]
                    k,n = len(countsVec),len(allcounts)
                    dfBetween,dfWithin = k-1,n-k
                    msr,mse = 0,0
                    for grp in countsVecAsArrays: msr += grp.size*(numpy.mean(grp)-grandmean)**2/float(dfBetween)
                    for grp,mn in zip(countsVecAsArrays,groupmeans): mse += numpy.sum((grp-mn)**2) 
                    mse /= float(dfWithin)
                    mse = mse+self.alpha ### moderation
                    Fmod = msr/float(mse)
                    Pmod = scipy.stats.f.sf(Fmod,dfBetween,dfWithin)
                    Fstat,pval = Fmod,Pmod
            pvals.append(pval)   
            Fstats.append(Fstat) 
            MSR.append(msr)
            MSE.append(mse)
            Rvs.append(Rv)

            # Update progress
            text = "Running Anova Method... %5.1f%%" % (100.0 * count / len(genes))
            self.progress_update(text, count)

        pvals = numpy.array(pvals)
        mask = numpy.isfinite(pvals)
        qvals = numpy.full(pvals.shape, numpy.nan)
        qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals[mask])[
            1
        ]  # BH, alpha=0.05

        msr, mse, f, p, q, statusMap = {},{},{},{},{},{}
        for i,rv in enumerate(Rvs):
          msr[rv],mse[rv],f[rv],p[rv],q[rv],statusMap[rv] = MSR[i],MSE[i],Fstats[i],pvals[i],qvals[i],status[i]
        return (msr, mse, f, p, q, statusMap)

    def calc_lf_cs(self, means, refs=[], pseudocount=1):
        if len(refs) == 0:
            refs = means  # if ref condition(s) not explicitly defined, use mean of all
        grandmean = numpy.mean(refs)
        lfcs = [math.log((x + pseudocount) / float(grandmean + pseudocount), 2) for x in means]
        return lfcs

    def Run(self):
        transit_tools.log("Starting Anova analysis")
        start_time = time.time()
        
        transit_tools.log("Getting Data")
        (sites, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(
            self.combined_wig
        )
        
        transit_tools.log("Normalizing using: %s" % self.normalization)
        (data, factors) = norm_tools.normalize_data(data, self.normalization)
        if self.winz:
            transit_tools.log("Winsorizing insertion counts")

        conditionsByFile, _, _, orderingMetadata = tnseq_tools.read_samples_metadata(
            self.metadata
        )
        conditions = self.wigs_to_conditions(conditionsByFile, filenamesInCombWig)

        conditionsList = self.select_conditions(
            conditions,
            self.included_conditions,
            self.excluded_conditions,
            orderingMetadata,
        )

        conditionNames = [conditionsByFile[f] for f in filenamesInCombWig]
        fileNames = filenamesInCombWig

        (
            data,
            fileNames,
            conditionNames,
            conditions,
            _,
            _,
        ) = self.filter_wigs_by_conditions3(  # in base.py
            data,
            fileNames,  # it looks like fileNames and conditionNames have to be parallel to data (vector of wigs)
            conditionNames,  # original Condition column in samples metadata file
            self.included_conditions,
            self.excluded_conditions,
            conditions=conditionNames,
        )  # this is kind of redundant for ANOVA, but it is here because condition, covars, and interactions could have been manipulated for ZINB

        genes = tnseq_tools.read_genes(self.annotation_path)

        TASiteindexMap = {TA: i for i, TA in enumerate(sites)}
        RvSiteindexesMap = tnseq_tools.rv_siteindexes_map(
            genes, TASiteindexMap, n_terminus=self.n_terminus, c_terminus=self.c_terminus
        )
        MeansByRv = self.means_by_rv(data, RvSiteindexesMap, genes, conditions)

        transit_tools.log("Running Anova")
        MSR, MSE, Fstats, pvals, qvals, run_status = self.run_anova(
            data, genes, MeansByRv, RvSiteindexesMap, conditions
        )

        transit_tools.log("Adding File: %s" % (self.output))
        file = open(self.output, "w")

        heads = (
            "Rv Gene TAs".split() +
            ["Mean_%s" % x for x in conditionsList] +
            ["LFC_%s" % x for x in conditionsList] +
            "MSR MSE+alpha Fstat Pval Padj".split() + 
            ["status"]
        )
        file.write("#Console: python3 %s\n" % " ".join(sys.argv))
        file.write("#parameters: normalization=%s, trimming=%s/%s%% (N/C), pseudocounts=%s, alpha=%s\n" % (self.normalization,self.n_terminus,self.c_terminus,self.pseudocount,self.alpha))
        file.write('#'+'\t'.join(heads)+EOL)
        for gene in genes:
            Rv = gene["rv"]
            if Rv in MeansByRv:
              means = [MeansByRv[Rv][c] for c in conditionsList]
              refs = [MeansByRv[Rv][c] for c in self.refs]
              LFCs = self.calc_lf_cs(means,refs,self.pseudocount)
              vals = ([Rv, gene["gene"], str(len(RvSiteindexesMap[Rv]))] +
                      ["%0.2f" % x for x in means] + 
                      ["%0.3f" % x for x in LFCs] + 
                      ["%f" % x for x in [MSR[Rv], MSE[Rv], Fstats[Rv], pvals[Rv], qvals[Rv]]] + [run_status[Rv]])
              file.write('\t'.join(vals)+EOL)
        file.close()
        transit_tools.log("Finished Anova analysis")
        transit_tools.log("Time: %0.1fs\n" % (time.time() - start_time))
    