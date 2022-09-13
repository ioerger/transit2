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

from pytransit.methods import analysis_base as base
import pytransit.tools.transit_tools as transit_tools
import pytransit.tools.tnseq_tools as tnseq_tools
import pytransit.tools.norm_tools as norm_tools
import pytransit.tools.stat_tools as stat_tools

# method_name = "griffin"


############# GUI ELEMENTS ##################

short_name = "griffin"
long_name = "Griffin"
short_desc = "Basic frequentist analysis of essentiality using gaps."
long_desc = "Analysis of gaps used in Griffin et al. 2011"
transposons = ["himar1"]
columns = [
    "Orf",
    "Name",
    "Desc",
    "k",
    "n",
    "r",
    "s",
    "t",
    "Expected Run",
    "p-value",
    "p-adjusted",
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
            GriffinMethod,
            GriffinGUI,
            [GriffinFile],
        )


################## FILE ###################


class GriffinFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#Griffin", columns)

    def get_header(self, path):
        ess = 0
        unc = 0
        non = 0
        short = 0
        with open(path) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                tmp = line.strip().split("\t")
                if float(tmp[-1]) < 0.05:
                    ess += 1
                else:
                    non += 1

        text = """Results:
    Essentials: %s
    Non-Essential: %s
            """ % (
            ess,
            non,
        )
        return text


################## GUI ###################


class GriffinGUI(base.AnalysisGUI):
    def define_panel(self, wxobj):
        self.wxobj = wxobj
        griffinPanel = wx.Panel(
            self.wxobj,
            wx.ID_ANY,
            wx.DefaultPosition,
            wx.DefaultSize,
            wx.TAB_TRAVERSAL,
        )

        griffinSection = wx.BoxSizer(wx.VERTICAL)

        griffinLabel = wx.StaticText(
            griffinPanel,
            wx.ID_ANY,
            u"griffin Options",
            wx.DefaultPosition,
            (120, -1),
            0,
        )
        griffinLabel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        griffinSection.Add(griffinLabel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        griffinSizer1 = wx.BoxSizer(wx.HORIZONTAL)

        griffinSection.Add(griffinSizer1, 1, wx.EXPAND, 5)

        griffinButton = wx.Button(
            griffinPanel,
            wx.ID_ANY,
            u"Run griffin",
            wx.DefaultPosition,
            wx.DefaultSize,
            0,
        )
        griffinSection.Add(griffinButton, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        griffinPanel.SetSizer(griffinSection)
        griffinPanel.Layout()
        griffinSection.Fit(griffinPanel)

        # Connect events
        griffinButton.Bind(wx.EVT_BUTTON, self.wxobj.RunMethod)

        self.panel = griffinPanel


########## CLASS #######################


class GriffinMethod(base.SingleConditionMethod):
    usage_string = """python3 %s griffin <comma-separated .wig files> <annotation .prot_table> <output file> [Optional Arguments]

        Optional Arguments:
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        -sC             :=  Include stop-codon (default is to ignore).
        -iN <float>     :=  Ignore TAs occuring at given fraction (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring at given fraction (as integer) of the C terminus. Default: -iC 0
        """ % sys.argv[0]
        
    def __init__(
        self,
        ctrldata,
        annotation_path,
        output_file,
        minread=1,
        replicates="Sum",
        normalization=None,
        LOESS=False,
        ignore_codon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
    ):

        base.SingleConditionMethod.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            ctrldata,
            annotation_path,
            output_file,
            replicates=replicates,
            normalization=normalization,
            LOESS=LOESS,
            ignore_codon=ignore_codon,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            wxobj=wxobj,
        )
        self.minread = minread

    @classmethod
    def from_gui(self, wxobj):
        """ """

        # Get Annotation file
        annotationPath = wxobj.annotation
        if not transit_tools.validate_annotation(annotationPath):
            return None

        # Get selected files
        ctrldata = wxobj.ctrlSelected()
        if not transit_tools.validate_control_datasets(ctrldata):
            return None

        # Validate transposon types
        if not transit_tools.validate_transposons_used(ctrldata, transposons):
            return None

        #
        minread = 1

        # Read the parameters from the wxPython widgets
        ignore_codon = True
        n_terminus = float(wxobj.globalNTerminusText.GetValue())
        c_terminus = float(wxobj.globalCTerminusText.GetValue())
        replicates = "Sum"
        normalization = None
        LOESS = False

        # Get output path
        name = transit_tools.basename(ctrldata[0])
        defaultFileName = "griffin_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path:
            return None
        output_file = open(output_path, "w")

        return self(
            ctrldata,
            annotationPath,
            output_file,
            minread,
            replicates,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
            wxobj,
        )

    @classmethod
    def from_args(self, args, kwargs):

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        minread = int(kwargs.get("m", 1))
        replicates = kwargs.get("r", "Sum")
        normalization = None
        LOESS = False
        ignore_codon = not kwargs.get("sC", False)
        n_terminus = float(kwargs.get("iN", 0.0))
        c_terminus = float(kwargs.get("iC", 0.0))

        return self(
            ctrldata,
            annotationPath,
            output_file,
            minread,
            replicates,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
        )

    def Run(self):

        transit_tools.log("Starting Griffin Method")
        start_time = time.time()

        # Get orf data
        transit_tools.log("Getting Data")

        (data, position) = transit_tools.get_validated_data(
            self.ctrldata, wxobj=self.wxobj
        )
        (K, N) = data.shape

        if self.normalization and self.normalization != "nonorm":
            transit_tools.log("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(
                data, self.normalization, self.ctrldata, self.annotation_path
            )

        G = tnseq_tools.Genes(
            self.ctrldata,
            self.annotation_path,
            minread=1,
            reps=self.replicates,
            ignore_codon=self.ignore_codon,
            n_terminus=self.n_terminus,
            c_terminus=self.c_terminus,
            data=data,
            position=position,
        )

        N = len(G)
        
        count = 0
        pins = G.global_theta()
        pnon = 1.0 - pins
        results = []
        for gene in G:
            if gene.n == 0:
                results.append([gene, 0.0, 1.000])
            else:
                B = 1.0 / math.log(1.0 / pnon)
                u = math.log(gene.n * pins, 1.0 / pnon)
                exprun = tnseq_tools.expected_runs(gene.n, pnon)
                pval = 1.0 - tnseq_tools.gumbel_cdf(gene.r, u, B)
                results.append([gene, exprun, pval])
            
            percentage = (100.0 * (count + 1) / (N))
            text = "Running Griffin Method... %5.1f%%" % percentage
            progress_update(text, percentage)
            count += 1

        pval = [row[-1] for row in results]
        padj = stat_tools.bh_fdr_correction(pval)
        for i in range(len(results)):
            results[i].append(padj[i])
        results.sort()

        self.output.write("#Griffin\n")
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
        self.output.write("#%s\n" % "\t".join(columns))

        for (gene, exprun, pval, padj) in results:
            self.output.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.1f\t%1.5f\t%1.5f\n"
                % (
                    gene.orf,
                    gene.name,
                    gene.desc,
                    gene.k,
                    gene.n,
                    gene.r,
                    gene.s,
                    gene.t,
                    exprun,
                    pval,
                    padj,
                )
            )

        self.output.close()

        transit_tools.log("")  # Printing empty line to flush stdout
        transit_tools.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)
        self.finish()
        transit_tools.log("Finished Griffin Method")

