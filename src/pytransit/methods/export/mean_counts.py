from pytransit.components.parameter_panel import panel, progress_update
import sys

from pytransit.tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub

import os
import time
import math
import random
import numpy
import scipy.stats
import datetime

from pytransit.methods import export_base as base
import pytransit
from pytransit.tools import transit_tools
from pytransit.tools import tnseq_tools
from pytransit.tools import norm_tools
from pytransit.tools import stat_tools


############# Description ##################

short_name = "mean_counts"
long_name = "Method to export datasets in 'Mean Gene Counts' format."
description = "A method to export and normalized datasets in 'Mean Gene Counts' format."
label = "to Mean Gene Counts"
transposons = ["himar1", "tn5"]

############# Analysis Method ##############


class Export(base.TransitExport):
    def __init__(self):
        base.TransitExport.__init__(
            self,
            short_name,
            long_name,
            description,
            label,
            transposons,
            MeanCountsMethod,
            MeanCountsGUI,
        )


################# GUI ##################


class MeanCountsGUI(base.ExportGUI):
    def __init__(self):
        base.ExportGUI.__init__(self)


########## METHOD #######################


class MeanCountsMethod(base.SingleConditionMethod):
    """   
    Mean Gene Counts
 
    """

    def __init__(
        self,
        combined_wig,
        ctrldata,
        annotation_path,
        output_file,
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
            description,
            label,
            ctrldata,
            annotation_path,
            output_file,
            normalization=normalization,
            LOESS=LOESS,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            wxobj=wxobj,
        )

        self.combined_wig = combined_wig  # boolean, interprete ctrldata as combined_wig file or comma-separated list of wig files?

    @classmethod
    def from_gui(self, wxobj):
        """ """

        # Get Annotation file
        from pytransit.interfaces import gui, cli
        annotation_path = universal.annotation_path
        if not transit_tools.validate_annotation(annotation_path):
            return None

        # Get selected files
        ctrldata = wxobj.ctrlSelected()
        if not transit_tools.validate_control_datasets(ctrldata):
            return None

        # Validate transposon types
        if not transit_tools.validate_transposons_used(ctrldata, transposons):
            return None

        # Choose normalization method
        normalization = wxobj.choose_normalization()

        LOESS = False
        ignore_codon = True
        n_terminus = 0.0
        c_terminus = 0.0

        # Get output path
        defaultFileName = "gene_mean_counts_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path:
            return None
        output_file = open(output_path, "w")

        # could add a checkbox for combined_wig
        combined_wig = False

        return self(
            combined_wig,
            ctrldata,
            annotation_path,
            output_file,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
            wxobj,
        )

    @classmethod
    def from_args(self, args, kwargs):
        print("ARGS=" + str(args))
        print("KWARGS=" + str(kwargs))

        combined_wig = kwargs.get("c", False)
        ctrldata = args[0].split(",")
        annotation_path = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        normalization = kwargs.get("n", "TTR")
        LOESS = False
        ignore_codon = True
        n_terminus = 0.0
        c_terminus = 0.0

        return self(
            combined_wig,
            ctrldata,
            annotation_path,
            output_file,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
        )

    def Run(self):

        logging.log("Starting Gene Mean Counts Export")
        start_time = time.time()

        # Get orf data
        logging.log("Getting Data")
        if self.combined_wig:
            (position, fulldata, datasets) = tnseq_tools.read_combined_wig(
                self.ctrldata[0]
            )
        else:
            (fulldata, position) = tnseq_tools.CombinedWig.gather_wig_data(self.ctrldata)
        (fulldata, factors) = norm_tools.normalize_data(
            fulldata, self.normalization, self.ctrldata, self.annotation_path
        )
        position = position.astype(int)

        hash = transit_tools.get_pos_hash(self.annotation_path)
        rv2info = transit_tools.get_gene_info(self.annotation_path)

        logging.log("Normalizing")
        self.output.write("#Summarized to Mean Gene Counts with TRANSIT.\n")
        if self.normalization != "nonorm":
            self.output.write("#Reads normalized using '%s'\n" % self.normalization)
            if type(factors[0]) == type(0.0):
                self.output.write(
                    "#Normalization Factors: %s\n"
                    % "\t".join(["%s" % f for f in factors.flatten()])
                )
            else:
                self.output.write(
                    "#Normalization Factors: %s\n"
                    % " ".join([",".join(["%s" % bx for bx in b]) for b in factors])
                )

        self.output.write("#Files:\n")
        names = datasets if self.combined_wig else self.ctrldata
        for f in names:
            self.output.write("#%s\n" % f)

        K, Nsites = fulldata.shape
        # Get Gene objects
        if self.combined_wig:
            G = tnseq_tools.Genes(
                self.ctrldata,
                self.annotation_path,
                norm=self.normalization,
                data=fulldata,
                position=position,
            )
        else:
            G = tnseq_tools.Genes(
                self.ctrldata, self.annotation_path, norm=self.normalization
            )
        N = len(G)
        
        if self.combined_wig:
            dataset_header = "\t".join(datasets)
        else:
            dataset_header = "\t".join(
                [transit_tools.fetch_name(D) for D in self.ctrldata]
            )
        self.output.write("#Orf\tName\tNumber of TA sites\t%s\n" % dataset_header)
        for i, gene in enumerate(G):
            if gene.n > 0:
                data_str = "\t".join(["%1.2f" % (M) for M in numpy.mean(gene.reads, 1)])
            else:
                data_str = "\t".join(["%1.2f" % (Z) for Z in numpy.zeros(K)])
            self.output.write(
                "%s\t%s\t%s\t%s\n" % (gene.orf, gene.name, gene.n, data_str)
            )

            # Update progress
            percentage = (100.0 * i / N)
            text = "Running Export Method... %5.1f%%" % percentage
            progress_update(text, percentage)
        self.output.close()

        logging.log("")  # Printing empty line to flush stdout
        self.finish()
        logging.log("Finished Export")

    #

    usage_string = (
            """python %s export mean_counts <comma-separated .wig files>|<combined_wig> <annotation .prot_table> <output file> [-c]\n note: append -c if inputing a combined_wig file\n"""
            % (sys.argv[0])
        )
