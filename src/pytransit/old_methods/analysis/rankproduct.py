from pytransit.components.parameter_panel import panel, progress_update
import pytransit.components.results_area as results_area
import sys

from pytransit.tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub


import os
import time
import math
import random
import numpy
import scipy.stats
import datetime


from pytransit.old_methods import analysis_base as base
from pytransit.tools import transit_tools
from pytransit.tools import tnseq_tools
from pytransit.tools import norm_tools
from pytransit.tools import stat_tools


############# GUI ELEMENTS ##################

short_name = "rankproduct"
long_name = "Rank Product"
short_desc = "Rank Product test for determining conditional essentiality."
long_desc = "Differential Comparison based on ranks"
transposons = ["himar1", "tn5"]
columns = [
    "Orf",
    "Name",
    "Desc",
    "Sites",
    "Mean Ctrl",
    "Mean Exp",
    "log2FC",
    "Obs RP",
    "Expected RP",
    "q-value",
]


############# Analysis Method ##############


class Analysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            transposons,
            RankProductMethod,
            RankProductGUI,
            [RankProductFile],
        )


################## FILE ###################


class RankProductFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#RankProduct", columns)


############# GUI ##################


class RankProductGUI(base.AnalysisGUI):
    def define_panel(self, wxobj):
        self.wxobj = wxobj
        rankproductPanel = wx.Panel(
            self.wxobj,
            wx.ID_ANY,
            wx.DefaultPosition,
            wx.DefaultSize,
            wx.TAB_TRAVERSAL,
        )

        rankproductSizer = wx.BoxSizer(wx.VERTICAL)

        rankproductLabel = wx.StaticText(
            rankproductPanel,
            wx.ID_ANY,
            "rankproduct Options",
            wx.DefaultPosition,
            (160, -1),
            0,
        )
        rankproductLabel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        rankproductSizer.Add(
            rankproductLabel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
        )

        rankproductLabelSizer = wx.BoxSizer(wx.VERTICAL)

        mainSizer1 = wx.BoxSizer(wx.VERTICAL)

        # SAMPLES
        (
            rankproductSampleLabel,
            self.wxobj.rankproductSampleText,
            sampleSizer,
        ) = self.defineTextBox(
            rankproductPanel,
            "Samples:",
            "10000",
            "Number of samples to take when estimating the theoretical rankproduct distribution. Larger samples give more accurate estimates at the cost of computation time.",
        )
        mainSizer1.Add(sampleSizer, 1, wx.EXPAND, 5)

        # NORMALIZATION
        # Norm
        rankproductNormChoiceChoices = [
            "TTR",
            "nzmean",
            "totreads",
            "zinfnb",
            "quantile",
            "betageom",
            "nonorm",
        ]
        (
            rankproductNormLabel,
            self.wxobj.rankproductNormChoice,
            normSizer,
        ) = self.defineChoiceBox(
            rankproductPanel,
            "Normalization:",
            rankproductNormChoiceChoices,
            "Choice of normalization method. The default choice, 'TTR', normalizes datasets to have the same expected count (while not being sensative to outliers). Read documentation for a description other methods. ",
        )
        mainSizer1.Add(normSizer, 1, wx.EXPAND, 5)

        rankproductSizer.Add(mainSizer1, 1, wx.EXPAND, 5)

        rankproductButton = wx.Button(
            rankproductPanel,
            wx.ID_ANY,
            "Run rankproduct",
            wx.DefaultPosition,
            wx.DefaultSize,
            0,
        )
        rankproductSizer.Add(
            rankproductButton, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
        )

        rankproductPanel.SetSizer(rankproductSizer)
        rankproductPanel.Layout()
        rankproductSizer.Fit(rankproductPanel)

        # Connect events
        rankproductButton.Bind(wx.EVT_BUTTON, self.wxobj.RunMethod)

        self.panel = rankproductPanel


########## CLASS #######################


class RankProductMethod(base.DualConditionMethod):
    usage_string = """python3 %s rankproduct <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 100
        -n <string>     :=  Normalization method. Default: -n TTR
        -h              :=  Output histogram of the permutations for each gene. Default: Turned Off.
        -a              :=  Perform adaptive rankproduct. Default: Turned Off.
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring at given fraction (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring at given fraction (as integer) of the C terminus. Default: -iC 0
        """ % sys.argv[0]

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
        replicates="Sum",
        LOESS=False,
        ignore_codon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
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
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            wxobj=wxobj,
        )

        self.samples = samples
        self.adaptive = adaptive
        self.doHistogram = doHistogram

    @classmethod
    def from_gui(self, wxobj):
        """ """

        # Get Annotation file
        from pytransit.universal_data import universal
        annotation_path = universal.annotation_path
        if not transit_tools.validate_annotation(annotation_path):
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
        ignore_codon = True
        samples = int(wxobj.rankproductSampleText.GetValue())
        normalization = wxobj.rankproductNormChoice.GetString(
            wxobj.rankproductNormChoice.GetCurrentSelection()
        )
        replicates = "Sum"
        adaptive = False
        doHistogram = False

        n_terminus = float(wxobj.globalNTerminusText.GetValue())
        c_terminus = float(wxobj.globalCTerminusText.GetValue())
        LOESS = False

        # Get output path
        defaultFileName = "rankproduct_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path:
            return None
        output_file = open(output_path, "w")

        return self(
            ctrldata,
            expdata,
            annotation_path,
            output_file,
            normalization,
            samples,
            adaptive,
            doHistogram,
            replicates,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
            wxobj,
        )

    @classmethod
    def from_args(self, args, kwargs):

        ctrldata = args[0].split(",")
        expdata = args[1].split(",")
        annotation_path = args[2]
        output_path = args[3]
        output_file = open(output_path, "w")

        normalization = kwargs.get("n", "TTR")
        samples = int(kwargs.get("s", 100))
        adaptive = kwargs.get("a", False)
        doHistogram = kwargs.get("h", False)
        replicates = kwargs.get("r", "Sum")

        LOESS = kwargs.get("l", False)
        ignore_codon = True
        n_terminus = float(kwargs.get("iN", 0.00))
        c_terminus = float(kwargs.get("iC", 0.00))

        return self(
            ctrldata,
            expdata,
            annotation_path,
            output_file,
            normalization,
            samples,
            adaptive,
            doHistogram,
            replicates,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
        )

    def Run(self):

        logging.log("Starting rankproduct Method")
        start_time = time.time()

        Kctrl = len(self.ctrldata)
        Kexp = len(self.expdata)
        # Get orf data
        logging.log("Getting Data")
        (data, position) = transit_tools.get_validated_data(
            self.ctrldata + self.expdata, wxobj=self.wxobj
        )
        if self.normalization != "none":
            logging.log("Normalizing using: %s" % self.normalization)

            (data, factors) = norm_tools.normalize_data(
                data,
                self.normalization,
                self.ctrldata + self.expdata,
                self.annotation_path,
            )

        Gctrl = tnseq_tools.Genes(
            self.ctrldata + self.expdata,
            self.annotation_path,
            ignore_codon=self.ignore_codon,
            n_terminus=self.n_terminus,
            c_terminus=self.c_terminus,
            data=data[:Kctrl, :],
            position=position,
        )

        Gexp = tnseq_tools.Genes(
            self.ctrldata + self.expdata,
            self.annotation_path,
            ignore_codon=self.ignore_codon,
            n_terminus=self.n_terminus,
            c_terminus=self.c_terminus,
            data=data[Kctrl:, :],
            position=position,
        )

        Ngenes = len(Gctrl)

        # Get the average counts for all the genes, in each replicate
        meanCtrl = numpy.zeros((Kctrl, Ngenes))
        meanExp = numpy.zeros((Kexp, Ngenes))

        for i in range(Ngenes):
            if numpy.any(Gctrl[i].reads):
                meanCtrl[:, i] = numpy.mean(Gctrl[i].reads, 1)
            else:
                meanCtrl[:, i] = numpy.zeros(Kctrl)
            #
            if numpy.any(Gexp[i].reads):
                meanExp[:, i] = numpy.mean(Gexp[i].reads, 1)
            else:
                meanExp[:, i] = numpy.zeros(Kexp)

        # Calculate a logFC2 between Experimental and Control
        # Then calculates it's rank, and observed rankProduct
        logFC2 = numpy.log2((meanExp + 0.0001) / (meanCtrl + 0.0001))
        rank = numpy.array([scipy.stats.rankdata(Lvec) for Lvec in logFC2])
        obsRP = numpy.power(numpy.prod(rank, 0), 1.0 / Kctrl)

        permutations = numpy.zeros((self.samples, Ngenes))
        tempranks = scipy.array([numpy.arange(1, Ngenes + 1) for rep in range(Kctrl)])
        for s in range(self.samples):
            rankperm = numpy.array([numpy.random.permutation(tr) for tr in tempranks])
            permutations[s] = numpy.power(numpy.prod(rankperm, 0), 1.0 / Kctrl)

        rankRP = numpy.argsort(obsRP) + 1

        # rankproduct
        data = []
        count = 0
        
        for i, gene in enumerate(Gctrl):
            count += 1

            meanctrl = numpy.mean(Gctrl[i].reads)
            meanexp = numpy.mean(Gexp[i].reads)
            log2fc = numpy.log2((meanexp + 0.0001) / (meanctrl + 0.0001))
            countbetter = numpy.sum(permutations <= obsRP[i])

            pval = countbetter / float(self.samples * Ngenes)
            e_val = countbetter / float(self.samples)
            q_paper = e_val / float(rankRP[i])

            data.append(
                [
                    gene.orf,
                    gene.name,
                    gene.desc,
                    gene.n,
                    meanctrl,
                    meanexp,
                    log2fc,
                    obsRP[i],
                    e_val,
                    q_paper,
                    pval,
                ]
            )

            # Update Progress
            percentage = (100.0 * count / Ngenes)
            text = "Running rankproduct Method... %5.1f%%" % percentage
            progress_update(text, percentage)

        #
        logging.log("")  # Printing empty line to flush stdout
        logging.log("Performing Benjamini-Hochberg Correction")
        data.sort()
        q_bh = stat_tools.bh_fdr_correction([row[-1] for row in data])

        self.output.write("#RankProduct\n")
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
                "#GUI with: ctrldata=%s, annotation=%s, output=%s\n"
                % (
                    ",".join(self.ctrldata).encode("utf-8"),
                    self.annotation_path.encode("utf-8"),
                    self.output.name.encode("utf-8"),
                )
            )
        else:
            self.output.write("#Console: python3 %s\n" % " ".join(sys.argv))

        self.output.write("#Data: %s\n" % (",".join(self.ctrldata).encode("utf-8")))
        self.output.write(
            "#Annotation path: %s\n" % self.annotation_path.encode("utf-8")
        )
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#%s\n" % (columns))

        for i, row in enumerate(data):
            (
                orf,
                name,
                desc,
                n,
                mean1,
                mean2,
                log2FCgene,
                obsRPgene,
                e_val,
                q_paper,
                pval,
            ) = row
            self.output.write(
                "%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.8f\t%1.1f\t%1.8f\n"
                % (
                    orf,
                    name,
                    desc,
                    n,
                    mean1,
                    mean2,
                    log2FCgene,
                    obsRPgene,
                    e_val,
                    q_paper,
                )
            )
        self.output.close()

        logging.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)
        self.finish()
        logging.log("Finished rankproduct Method")

