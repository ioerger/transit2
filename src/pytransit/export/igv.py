import sys

from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub


import os
import time
import math
import random
import numpy
import scipy.stats
import datetime

from pytransit.export import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools


############# Description ##################

short_name = "igv"
long_name = "Method to export datasets in 'IGV' format."
description = "A method to export and normalized datasets in 'IGV' format."
label = "to IGV"
transposons = ["himar1", "tn5"]

############# Analysis Method ##############


class IGVExport(base.TransitExport):
    def __init__(self):
        base.TransitExport.__init__(
            self,
            short_name,
            long_name,
            description,
            label,
            transposons,
            IGVMethod,
            IGVGUI,
        )


################# GUI ##################


class IGVGUI(base.ExportGUI):
    def __init__(self):
        base.ExportGUI.__init__(self)


########## METHOD #######################


class IGVMethod(base.SingleConditionMethod):
    """   
    IGV
 
    """

    def __init__(
        self,
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

        # Choose normalization method
        normalization = wxobj.chooseNormalization()

        LOESS = False
        ignore_codon = True
        n_terminus = 0.0
        c_terminus = 0.0

        # Get output path
        defaultFileName = "igv_output.txt"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path:
            return None
        output_file = open(output_path, "w")

        return self(
            ctrldata,
            annotationPath,
            output_file,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
            wxobj,
        )

    @classmethod
    def from_args(self, rawargs):
        (args, kwargs) = transit_tools.clean_args(rawargs)

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        normalization = kwargs.get("n", "TTR")
        LOESS = False
        ignore_codon = True
        n_terminus = 0.0
        c_terminus = 0.0

        return self(
            ctrldata,
            annotationPath,
            output_file,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
        )

    def Run(self):

        transit_tools.log("Starting IGV Export")
        start_time = time.time()

        # Get orf data
        transit_tools.log("Getting Data")
        (fulldata, position) = tnseq_tools.CombinedWig.gather_wig_data(self.ctrldata)
        (fulldata, factors) = norm_tools.normalize_data(
            fulldata, self.normalization, self.ctrldata, self.annotation_path
        )
        position = position.astype(int)

        hash = transit_tools.get_pos_hash(self.annotation_path)
        rv2info = transit_tools.get_gene_info(self.annotation_path)

        transit_tools.log("Normalizing")
        self.output.write("#Converted to IGV with TRANSIT.\n")
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
        for f in self.ctrldata:
            self.output.write("#%s\n" % f)

        dataset_str = "\t".join([transit_tools.fetch_name(F) for F in self.ctrldata])
        self.output.write("#Chromosome\tStart\tEnd\tFeature\t%s\tTAs\n" % dataset_str)
        chrom = transit_tools.fetch_name(self.annotation_path)

        (K, N) = fulldata.shape
        
        for i, pos in enumerate(position):
            self.output.write(
                "%s\t%s\t%s\tTA%s\t%s\t1\n"
                % (
                    chrom,
                    position[i],
                    position[i] + 1,
                    position[i],
                    "\t".join(["%1.1f" % fulldata[j][i] for j in range(len(fulldata))]),
                )
            )

            # Update progress
            text = "Running Export Method... %5.1f%%" % (100.0 * i / N)
            self.progress_update(text, i)
        self.output.close()

        transit_tools.log("")  # Printing empty line to flush stdout
        self.finish()
        transit_tools.log("Finished Export")

    #

    @classmethod
    def usage_string(self):
        return (
            """python %s export igv <comma-separated .wig files> <annotation .prot_table> <output file>"""
            % (sys.argv[0])
        )


if __name__ == "__main__":

    (args, kwargs) = transit_tools.clean_args(sys.argv[1:])

    print("ARGS:", args)
    print("KWARGS:", kwargs)

    G = Example.from_args(sys.argv[1:])

    print(G)
    G.console_message("Printing the member variables:")
    G.print_members()

    print("")
    print("Running:")

    G.Run()
