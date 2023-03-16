from pytransit.components.parameter_panel import panel, progress_update
import pytransit.components.results_area as results_area
import sys

from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub


import os
import time
import math
import random
import numpy
import scipy.stats
import datetime

from pytransit.old_methods import analysis_base as base
from pytransit.specific_tools import transit_tools
from pytransit.specific_tools import tnseq_tools
from pytransit.specific_tools import norm_tools
from pytransit.specific_tools import stat_tools
from pytransit.globals import logging

# method_name = "example"


############# GUI ELEMENTS ##################

short_name = "tn5gaps"
long_name = "Tn5 Gaps"
short_desc = "Analysis of essentiality on gaps in entire genome (Tn5)."
long_desc = "A analysis method based on the extreme value (Gumbel) distribution that considers longest runs over the whole genome instead of individual genes."
transposons = ["tn5"]
columns = [
    "Orf",
    "Name",
    "Desc",
    "k",
    "n",
    "r",
    "ovr",
    "lenovr",
    "pval",
    "padj",
    "call",
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
            Tn5GapsMethod,
            Tn5GapsGUI,
            [Tn5GapsFile],
        )


################## FILE ###################


class Tn5GapsFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#Tn5 Gaps", columns)

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
                if tmp[-1] == "Essential":
                    ess += 1
                if tmp[-1] == "Non-essential":
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


class Tn5GapsGUI(base.AnalysisGUI):
    def define_panel(self, wxobj):
        self.wxobj = wxobj
        tn5GapsPanel = wx.Panel(
            self.wxobj,
            wx.ID_ANY,
            wx.DefaultPosition,
            wx.DefaultSize,
            wx.TAB_TRAVERSAL,
        )

        tn5GapsSection = wx.BoxSizer(wx.VERTICAL)

        tn5GapsLabel = wx.StaticText(
            tn5GapsPanel,
            wx.ID_ANY,
            "Tn5 Gaps Options",
            wx.DefaultPosition,
            (150, -1),
            0,
        )
        tn5GapsLabel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        tn5GapsSection.Add(tn5GapsLabel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        mainSizer1 = wx.BoxSizer(wx.VERTICAL)

        # Min Read
        tn5GapsReadChoiceChoices = ["1", "2", "3", "4", "5"]
        (
            tn5GapsReadLabel,
            self.wxobj.tn5GapsReadChoice,
            readSizer,
        ) = self.defineChoiceBox(
            tn5GapsPanel,
            "Minimum Read:",
            tn5GapsReadChoiceChoices,
            "This is the minimum number of reads to consider a 'true' insertion. Value of 1 will consider all insertions. Larger values allow the method to ignore spurious insertions which might interrupt a run of non-insertions. Noisy datasets or those with many replicates can beneffit from increasing this.",
        )
        mainSizer1.Add(readSizer, 1, wx.EXPAND, 5)

        # Replicates
        tn5GapsRepChoiceChoices = ["Sum", "Mean"]
        (tn5GapsRepLabel, self.wxobj.tn5GapsRepChoice, repSizer) = self.defineChoiceBox(
            tn5GapsPanel,
            "Replicates:",
            tn5GapsRepChoiceChoices,
            "Determines how to handle replicates, and their read-counts. When using many replicates, summing read-counts may make spurious counts appear to be significantly large and interrupt a run of non-insertions.",
        )
        mainSizer1.Add(repSizer, 1, wx.EXPAND, 5)

        tn5GapsSection.Add(mainSizer1, 1, wx.EXPAND, 5)

        tn5GapsButton = wx.Button(
            tn5GapsPanel,
            wx.ID_ANY,
            "Run Tn5Gaps",
            wx.DefaultPosition,
            wx.DefaultSize,
            0,
        )
        tn5GapsSection.Add(tn5GapsButton, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        tn5GapsPanel.SetSizer(tn5GapsSection)
        tn5GapsPanel.Layout()
        tn5GapsSection.Fit(tn5GapsPanel)

        # Connect events
        tn5GapsButton.Bind(wx.EVT_BUTTON, wxobj.RunMethod)

        self.panel = tn5GapsPanel


########## CLASS #######################


class Tn5GapsMethod(base.SingleConditionMethod):
    """   
    Example
 
    """

    def __init__(
        self,
        ctrldata,
        annotation_path,
        output_file,
        replicates="Sum",
        normalization=None,
        LOESS=False,
        ignore_codon=True,
        minread=1,
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
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            wxobj=wxobj,
        )
        self.minread = minread

    @classmethod
    def from_gui(self, wxobj):
        """ """
        # Get Annotation file
        from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
        annotation_path = gui.annotation_path
        if not transit_tools.validate_annotation(annotation_path):
            return None

        # Get selected files
        ctrldata = wxobj.ctrlSelected()
        if not transit_tools.validate_control_datasets(ctrldata):
            return None

        # Validate transposon types
        types = tnseq_tools.get_file_types(ctrldata)
        if "himar1" in types:
            answer = logging.warn(
                "Warning: One of the selected wig files looks like a Himar1 dataset. This method is designed to work on Tn5 wig files. Proceeding will fill in missing data with zeroes. Click OK to continue."
            )
            if answer == wx.ID_CANCEL:
                return None

        # Read the parameters from the wxPython widgets
        ignore_codon = True
        minread = int(
            wxobj.tn5GapsReadChoice.GetString(
                wxobj.tn5GapsReadChoice.GetCurrentSelection()
            )
        )
        n_terminus = float(wxobj.globalNTerminusText.GetValue())
        c_terminus = float(wxobj.globalCTerminusText.GetValue())
        replicates = wxobj.tn5GapsRepChoice.GetString(
            wxobj.tn5GapsRepChoice.GetCurrentSelection()
        )
        normalization = None
        LOESS = False

        # Get output path
        name = transit_tools.basename(ctrldata[0])
        defaultFileName = "tn5_gaps_output_m%d_r%s.dat" % (minread, replicates)

        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path:
            return None
        output_file = open(output_path, "w")

        return self(
            ctrldata,
            annotation_path,
            output_file,
            replicates,
            normalization,
            LOESS,
            ignore_codon,
            minread,
            n_terminus,
            c_terminus,
            wxobj,
        )

    @classmethod
    def from_args(self, args, kwargs):

        ctrldata = args[0].split(",")
        annotation_path = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        replicates = kwargs.get("r", "Sum")
        minread = int(kwargs.get("m", 1))
        normalization = None
        LOESS = False
        ignore_codon = True
        n_terminus = float(kwargs.get("iN", "0"))
        c_terminus = float(kwargs.get("iC", "0"))

        return self(
            ctrldata,
            annotation_path,
            output_file,
            replicates,
            normalization,
            LOESS,
            ignore_codon,
            minread,
            n_terminus,
            c_terminus,
        )

    def Run(self):
        logging.log("Starting Tn5 gaps method")
        start_time = time.time()

        logging.log("Loading data (May take a while)")

        # Combine all wigs
        (data, position) = transit_tools.get_validated_data(
            self.ctrldata, wxobj=self.wxobj
        )
        combined = tnseq_tools.combine_replicates(data, method=self.replicates)
        combined[combined < self.minread] = 0
        counts = combined
        counts[counts > 0] = 1
        num_sites = counts.size

        genes_obj = tnseq_tools.Genes(
            self.ctrldata,
            self.annotation_path,
            ignore_codon=self.ignore_codon,
            n_terminus=self.n_terminus,
            c_terminus=self.c_terminus,
            data=data,
            position=position,
        )

        pins = numpy.mean(counts)
        pnon = 1.0 - pins

        # Calculate stats of runs
        exprunmax = tnseq_tools.expected_runs(num_sites, pnon)
        varrun = tnseq_tools.variance_run(num_sites, pnon)
        stddevrun = math.sqrt(varrun)
        exp_cutoff = exprunmax + 2 * stddevrun

        # Get the runs
        logging.log("Identifying non-insertion runs in genome")
        run_arr = tnseq_tools.runs_w_info(counts)
        pos_hash = transit_tools.get_pos_hash(self.annotation_path)

        # Finally, calculate the results
        logging.log("Running Tn5 gaps method")
        results_per_gene = {}
        for gene in genes_obj.genes:
            results_per_gene[gene.orf] = [
                gene.orf,
                gene.name,
                gene.desc,
                gene.k,
                gene.n,
                gene.r,
                0,
                0,
                1,
            ]

        N = len(run_arr)
        count = 0
        accum = 0
        
        for run in run_arr:
            accum += run["length"]
            count += 1
            genes = tnseq_tools.get_genes_in_range(pos_hash, run["start"], run["end"])
            for gene_orf in genes:
                gene = genes_obj[gene_orf]  # bug fix: moved this up
                start, end = gene.start, gene.end
                a, b = self.n_terminus, self.c_terminus
                if gene.strand == "-":
                    a, b = b, a
                start = start + int((end - start) * (a / 100.0))
                end = end - int((end - start) * (b / 100.0))

                inter_sz = (
                    self.intersect_size([run["start"], run["end"]], [start, end]) + 1
                )
                percent_overlap = self.calc_overlap(
                    [run["start"], run["end"]], [start, end]
                )
                run_len = run["length"]
                B = 1.0 / math.log(1.0 / pnon)
                u = math.log(num_sites * pins, 1.0 / pnon)
                pval = 1.0 - tnseq_tools.gumbel_cdf(run["length"], u, B)

                curr_val = results_per_gene[gene.orf]
                curr_inter_sz = curr_val[6]
                curr_len = curr_val[7]
                if inter_sz > curr_inter_sz:
                    results_per_gene[gene.orf] = [
                        gene.orf,
                        gene.name,
                        gene.desc,
                        gene.k,
                        gene.n,
                        gene.r,
                        inter_sz,
                        run_len,
                        pval,
                    ]

            # Update Progress
            if count % 10000 == 0:
                percent = (100.0 * count / N)
                text = "Running Tn5Gaps method... %1.1f%%" % percent
                progress_update(text, percent)

        data = list(results_per_gene.values())
        exp_run_len = float(accum) / N

        min_sig_len = float("inf")
        sig_genes_count = 0
        pval = [row[-1] for row in data]
        padj = stat_tools.bh_fdr_correction(pval)
        for i in range(len(data)):
            if padj[i] < 0.05:
                sig_genes_count += 1
                min_sig_len = min(min_sig_len, data[i][-2])
            data[i].append(padj[i])
            data[i].append("Essential" if padj[i] < 0.05 else "Non-essential")
            # (data[i][0], data[i][1], data[i][2], data[i][3], data[i][4], data[i][5], data[i][6], data[i][7], data[i][8], padj[i], 'Essential' if padj[i] < 0.05 else 'Non-essential')
        data.sort(key=lambda l: l[0])

        # Output results
        self.output.write("#Tn5 Gaps\n")
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
        self.output.write("#Essential gene count: %d\n" % (sig_genes_count))
        self.output.write("#Minimum reads: %d\n" % (self.minread))
        self.output.write("#Replicate combination method: %s\n" % (self.replicates))
        self.output.write("#Insertion density: %0.3f\n" % (pins))
        self.output.write("#Mean run length: %0.1f\n" % (exp_run_len))
        self.output.write("#Expected max run length: %0.1f\n" % (exprunmax))
        self.output.write("#Minimum significant run length: %d\n" % (min_sig_len))
        self.output.write("#%s\n" % "\t".join(columns))
        # self.output.write("#Orf\tName\tDesc\tk\tn\tr\tovr\tlenovr\tpval\tpadj\tcall\n")

        for res in data:
            self.output.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%1.5f\t%1.5f\t%s\n"
                % (
                    res[0],
                    res[1],
                    res[2],
                    res[3],
                    res[4],
                    res[5],
                    res[6],
                    res[7],
                    res[8],
                    res[9],
                    res[10],
                )
            )
        self.output.close()

        logging.log("")  # Printing empty line to flush stdout
        logging.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)
        self.finish()
        logging.log("Finished Tn5Gaps Method")

    usage_string = """python3 %s tn5gaps <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    
        Optional Arguments:
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        -iN <float>     :=  Ignore TAs occurring within given percentage (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occurring within given percentage (as integer) of the C terminus. Default: -iC 0
        """ % sys.argv[0]

    def intersect_size(self, intv1, intv2):
        first = intv1 if intv1[0] < intv2[0] else intv2
        second = intv1 if first == intv2 else intv2

        if first[1] < second[0]:
            return 0

        right_ovr = min(first[1], second[1])
        left_ovr = max(first[0], second[0])

        return right_ovr - left_ovr

    def calc_overlap(self, run_interv, gene_interv):
        first = run_interv if run_interv[0] < gene_interv[0] else gene_interv
        second = run_interv if first == gene_interv else gene_interv

        if second[0] > first[0] and second[1] < first[1]:
            return 1.0

        intersect = self.intersect_size(run_interv, gene_interv)
        return float(intersect) / (gene_interv[1] - gene_interv[0])

