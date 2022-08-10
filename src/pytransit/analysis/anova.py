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

from pytransit.transit_tools import wx, pub, basename
from pytransit.analysis import base
import pytransit
import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools
import pytransit.components.results_area as results_area
from pytransit.core_data import universal
from pytransit.components.parameter_panel import panel as parameter_panel
from pytransit.components.parameter_panel import panel, progress_update
from pytransit.components.panel_helpers import make_panel, create_run_button, create_normalization_input, create_reference_condition_input, create_include_condition_list_input, create_exclude_condition_list_input, create_n_terminus_input, create_c_terminus_input, create_pseudocount_input, create_winsorize_input, create_alpha_input, create_button

############# GUI ELEMENTS ##################

class Analysis:
    identifier  = "#Anova"
    short_name  = "anova"
    long_name   = "ANOVA"
    short_desc  = "Perform Anova analysis"
    long_desc   = """Perform Anova analysis"""
    transposons = [ "himar1", "tn5" ]
    
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
    
    gui=None
    method=None
    
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
    )
    
    def __init__(self):
        self.full_name        = f"[{self.short_name}]  -  {self.short_desc}"
        self.transposons_text = transit_tools.get_transposons_text(self.transposons)
        self.method           = Method
        self.gui              = GUI()
        self.filetypes        = [File]
    
    def __str__(self):
        return f"""
            Analysis Method:
                Short Name:  {self.short_name}
                Long Name:   {self.long_name}
                Short Desc:  {self.short_desc}
                Long Desc:   {self.long_desc}
                Method:      {self.method}
                GUI:         {self.gui}
        """.replace('\n            ','\n').strip()

@transit_tools.ResultsFile
class File(Analysis):
    def __init__(self, path):
        self.wxobj = None
        self.path  = path
        self.table_values = dict(
            name=basename(self.path),
            type=Analysis.identifier,
            path=self.path,
            __define_panel=self.define_panel, # part of the row, but isn't shown on the table. Name is used by the results area
        )
    
    @staticmethod
    def can_load(path):
        row = 0
        data = []
        shown_error = False
        has_correct_identifier = False
        with open(path) as in_file:
            for line in in_file:
                if line.startswith("#"):
                    if line.startswith(Analysis.identifier):
                        has_correct_identifier = True
                    continue
                tmp = line.split("\t")
                tmp[-1] = tmp[-1].strip()
                try:
                    rowdict = dict([(Analysis.columns[i], tmp[i]) for i in range(len(Analysis.columns))])
                except Exception as e:
                    print(e)
                    return False
                data.append((row, rowdict))
                row += 1
        return has_correct_identifier
    
    def define_panel(self):
        self.panel = make_panel()
        
        self.value_getters = LazyDict()
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        if True:
            # 
            # heatmap options
            # 
            @create_button(self.panel, main_sizer, label="Generate Heatmap")
            def _(event):
                print(f'''event = {event}''')
        
        # 
        # finish panel setup
        # 
        parameter_panel.set_panel(self.panel)
        self.panel.SetSizer(main_sizer)
        self.panel.Layout()
        main_sizer.Fit(self.panel)
        
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

    def create_heatmap(filetype, infile, outfile, topk, qval, low_mean_filter):
        headers = None
        data, hits = [], []
        n = -1  # number of conditions

        for line in open(infile):
            w = line.rstrip().split("\t")
            if line[0] == "#" or (
                "pval" in line and "padj" in line
            ):  # check for 'pval' for backwards compatibility
                headers = w
                continue  # keep last comment line as headers
            # assume first non-comment line is header
            if n == -1:
                # ANOVA header line has names of conditions, organized as 3+2*n+3 (2 groups (means, LFCs) X n conditions)
                # ZINB header line has names of conditions, organized as 3+4*n+3 (4 groups X n conditions)
                if filetype == "anova":
                    n = int((len(w) - 6) / 2)
                elif filetype == "zinb":
                    n = int((len(headers) - 6) / 4)
                headers = headers[3 : 3 + n]
                headers = [x.replace("Mean_", "") for x in headers]
            else:
                means = [
                    float(x) for x in w[3 : 3 + n]
                ]  # take just the columns of means
                lfcs = [
                    float(x) for x in w[3 + n : 3 + n + n]
                ]  # take just the columns of LFCs
                qval = float(w[-2])
                data.append((w, means, lfcs, qval))

        data.sort(key=lambda x: x[-1])
        hits, LFCs = [], []
        for k, (w, means, lfcs, qval) in enumerate(data):
            if (topk == -1 and qval < qval) or (
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
        
        
        heatmapFunc = make_heatmapFunc()
        heatmapFunc(df, StrVector(genenames), outfile)
    

class GUI:
    def __init__(self, *args, **kwargs):
        Analysis.gui = self
        self.wxobj = None
        self.panel = None
    
    def define_panel(self, _):
        self.panel = make_panel()
        
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
            self.value_getters.included_conditions    = create_include_condition_list_input(self.panel, main_sizer)
            self.value_getters.excluded_conditions    = create_exclude_condition_list_input(self.panel, main_sizer)
            self.value_getters.reference_condition    = create_reference_condition_input(self.panel, main_sizer)
            self.value_getters.n_terminus             = create_n_terminus_input(self.panel, main_sizer)
            self.value_getters.c_terminus             = create_c_terminus_input(self.panel, main_sizer)
            self.value_getters.normalization          = create_normalization_input(self.panel, main_sizer)
            self.value_getters.pseudocount            = create_pseudocount_input(self.panel, main_sizer)
            self.value_getters.alpha                  = create_alpha_input(self.panel, main_sizer)
            self.value_getters.winz                   = create_winsorize_input(self.panel, main_sizer)
            self.value_getters.refs                   = lambda *args: [] if self.value_getters.reference_condition() == "[None]" else [ self.value_getters.reference_condition() ]
            
            create_run_button(self.panel, main_sizer)
            
        parameter_panel.set_panel(self.panel)
        self.panel.SetSizer(main_sizer)
        self.panel.Layout()
        main_sizer.Fit(self.panel)


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
        Analysis.method = self
        
        self.pseudocount         = None
        self.alpha               = None
        self.refs                = None
        self.winz                = None
        self.normalization       = None
        self.ctrldata            = None
        self.expdata             = None
        self.annotation_path     = None
        self.LOESS               = None
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
            Analysis.short_name,
            Analysis.long_name,
            Analysis.short_desc,
            Analysis.long_desc,
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
            Analysis.inputs.combined_wig = wig_group.cwig.path
            Analysis.inputs.metadata     = wig_group.metadata.path
            
            # 
            # get annotation
            # 
            Analysis.inputs.annotation = universal.session_data.annotation
            # FIXME: enable this once I get a valid annotation file example
            # if not transit_tools.validate_annotation(Analysis.inputs.annotation):
            #     return None
            
            # 
            # setup custom inputs
            # 
            for each_key, each_getter in Analysis.gui.value_getters.items():
                try:
                    Analysis.inputs[each_key] = each_getter()
                except Exception as error:
                    raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
            
            # 
            # save result files
            # 
            Analysis.inputs.output_file = gui_tools.ask_for_output_file_path(
                default_file_name="test1_output.dat",
                output_extensions=u'Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
            )
            if not Analysis.inputs.output_file:
                return None

            return self(*Analysis.inputs.values())

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
                raise Exception(f'''flag unrecognized: {arg}''')
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
            raise Exception(f'''Only 2 conditions expected, but got: {included_conditions}''')
            sys.exit(0)
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

    def calculate_anova(self, data, genes, MeansByRv, RvSiteindexesMap, conditions):
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
            percentage = 100.0 * count / len(genes)
            progress_update(f"Running Anova Method... {percentage:5.1f}%", percentage)

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
        
        # 
        # get data
        # 
        transit_tools.log("Getting Data")
        if True:
            sites, data, filenamesInCombWig = tnseq_tools.read_combined_wig(self.combined_wig)
            
            transit_tools.log(f"Normalizing using: {self.normalization}")
            data, factors = norm_tools.normalize_data(data, self.normalization)

            if self.winz: transit_tools.log("Winsorizing insertion counts")
            conditionsByFile, _, _, orderingMetadata = tnseq_tools.read_samples_metadata(self.metadata)
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
        
        # 
        # process data
        # 
        if True:
            TASiteindexMap = {TA: i for i, TA in enumerate(sites)}
            RvSiteindexesMap = tnseq_tools.rv_siteindexes_map(
                genes, TASiteindexMap, n_terminus=self.n_terminus, c_terminus=self.c_terminus
            )
            MeansByRv = self.means_by_rv(data, RvSiteindexesMap, genes, conditions)

            transit_tools.log("Running Anova")
            MSR, MSE, Fstats, pvals, qvals, run_status = self.calculate_anova(
                data, genes, MeansByRv, RvSiteindexesMap, conditions
            )
        
        # 
        # write output
        # 
        transit_tools.log(f"Adding File: {self.output}")
        results_area.add(self.output)
        if True:
            file = open(self.output, "w")

            heads = (
                "Rv Gene TAs".split() +
                ["Mean_%s" % x for x in conditionsList] +
                ["LFC_%s" % x for x in conditionsList] +
                "MSR MSE+alpha Fstat Pval Padj".split() + 
                ["status"]
            )
            file.write(Analysis.identifier+"\n")
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
    

def create_heatmap(filetype, infile, outfile, topk, qval, low_mean_filter):
    headers = None
    data, hits = [], []
    n = -1  # number of conditions

    for line in open(infile):
        w = line.rstrip().split("\t")
        if line[0] == "#" or (
            "pval" in line and "padj" in line
        ):  # check for 'pval' for backwards compatibility
            headers = w
            continue  # keep last comment line as headers
        # assume first non-comment line is header
        if n == -1:
            # ANOVA header line has names of conditions, organized as 3+2*n+3 (2 groups (means, LFCs) X n conditions)
            # ZINB header line has names of conditions, organized as 3+4*n+3 (4 groups X n conditions)
            if filetype == "anova":
                n = int((len(w) - 6) / 2)
            elif filetype == "zinb":
                n = int((len(headers) - 6) / 4)
            headers = headers[3 : 3 + n]
            headers = [x.replace("Mean_", "") for x in headers]
        else:
            means = [
                float(x) for x in w[3 : 3 + n]
            ]  # take just the columns of means
            lfcs = [
                float(x) for x in w[3 + n : 3 + n + n]
            ]  # take just the columns of LFCs
            qval = float(w[-2])
            data.append((w, means, lfcs, qval))

    data.sort(key=lambda x: x[-1])
    hits, LFCs = [], []
    for k, (w, means, lfcs, qval) in enumerate(data):
        if (topk == -1 and qval < qval) or (
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
    
    
    heatmapFunc = make_heatmapFunc()
    heatmapFunc(df, StrVector(genenames), outfile)