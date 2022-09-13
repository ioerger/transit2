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


############# Description ##################

short_name = "example"
long_name = "Example"
short_desc = "Example method that calculates mean read-counts per gene."
long_desc = "A method made to serve as an example to implementing other methods."
transposons = ["himar1", "tn5"]
columns = ["Orf", "Name", "Desc", "k", "n", "mean", "nzmean"]

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
            Method,
            GUI,
            [File],
        )


################## FILE ###################


class File(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#Example", columns)

    def get_header(self, path):
        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
        return text


################# GUI ##################


class GUI(base.AnalysisGUI):
    def __init__(self):
        base.AnalysisGUI.__init__(self)


########## METHOD #######################


class Method(base.SingleConditionMethod):
    usage_string = """python3 %s example <comma-separated .wig files> <annotation .prot_table> <output file>""" % (sys.argv[0])

    def __init__(
        self,
        ctrldata,
        annotation_path,
        output_file,
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
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            wxobj=wxobj,
        )

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

        # Read the parameters from the wxPython widgets
        ignore_codon = True
        n_terminus = float(wxobj.globalNTerminusText.GetValue())
        c_terminus = float(wxobj.globalCTerminusText.GetValue())
        replicates = "Sum"
        normalization = None
        LOESS = False

        # Get output path
        defaultFileName = "example_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path:
            return None
        output_file = open(output_path, "w")

        return self(
            ctrldata,
            annotationPath,
            output_file,
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

        replicates = "Sum"
        normalization = None
        LOESS = False
        ignore_codon = True
        n_terminus = 0.0
        c_terminus = 0.0

        return self(
            ctrldata,
            annotationPath,
            output_file,
            replicates,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
        )

    def Run(self):

        transit_tools.log("Starting Example Method")
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

        data = []
        N = len(G)
        count = 0
        
        for gene in G:
            count += 1
            if gene.n == 0:
                mean = 0.0
            else:
                mean = numpy.mean(gene.reads)

            if gene.k == 0:
                nzmean = 0.0
            else:
                nzmean = numpy.sum(gene.reads) / float(gene.k)

            data.append(
                "%s\t%s\t%s\t%s\t%s\t%1.2f\t%1.2f\n"
                % (gene.orf, gene.name, gene.desc, gene.k, gene.n, mean, nzmean)
            )

            # Update Progress
            percent = (100.0 * count / N)
            text = "Running Example Method... %5.1f%%" % percent
            progress_update(text, percent)

        self.output.write("#Example\n")
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

        data.sort()
        for line in data:
            self.output.write(line)
        self.output.close()

        transit_tools.log("")  # Printing empty line to flush stdout
        transit_tools.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)
        self.finish()
        transit_tools.log("Finished Example Method")

    

