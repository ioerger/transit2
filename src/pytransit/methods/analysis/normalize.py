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

short_name = "normalize"
long_name = "Normalize"
short_desc = "Normalization method"
long_desc = "Method for normalizing datasets."
transposons = ["himar1", "tn5"]
columns = ["Position", "Reads", "Genes"]


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
            NormalizeMethod,
            NormalizeGUI,
            [NormalizeFile],
        )


################## FILE ###################

# maybe this is only needed if it is going to get loaded into the output-files panel - TRI


class NormalizeFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#CombinedWig", columns)

    def get_header(self, path):
        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
        return text


################# GUI ##################


class NormalizeGUI(base.AnalysisGUI):
    def __init__(self):
        base.AnalysisGUI.__init__(self)


########## METHOD #######################


class NormalizeMethod(base.SingleConditionMethod):
    usage_string = """
        python3 %s normalize <input.wig> <output.wig> [-n TTR|betageom]
        ---
        OR
        ---
        python3 %s normalize -c <input combined_wig> <output.wig> [-n TTR|betageom]

            Optional Arguments:
            -n <string>     :=  Normalization method. Default: -n TTR
        """ % (sys.argv[0], sys.argv[0],)

    def __init__(self, infile, outfile, normalization):
        ctrldata = [infile]
        annotation_path = ""
        output_file = outfile
        replicates = "Sum"
        normalization = normalization
        LOESS = False
        ignore_codon = True
        n_terminus = 0.0
        c_terminus = 0.0
        wxobj = None
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
    def from_args(self, args, kwargs):
        from pytransit.tools.console_tools import InvalidArgumentException
        
        isCombinedWig = "c" in kwargs
        if (not isCombinedWig and len(args) < 2) or (isCombinedWig and len(args) < 1):
            raise InvalidArgumentException("Must provide all necessary arguments")
        if isCombinedWig:
            self.infile = kwargs.get("c")  # only 1 input wig file
            self.outfile = args[0]  # if no arg give, could print to screen
        else:
            self.infile = args[0]  # only 1 input wig file
            self.outfile = args[1]  # if no arg give, could print to screen
        self.normalization = kwargs.get(
            "n", "TTR"
        )  # check if it is a legal method name
        self.combined_wig = isCombinedWig

        return self(self.infile, self.outfile, self.normalization)

    def Run(self):

        logging.log("Starting Normalization")
        start_time = time.time()

        infile = self.infile
        outputPath = (
            self.outfile
        )  # output file exists, should I require -overwrite flag?

        # determine ref genome from first; assume they are all the same; assume wigs have 2 header lines
        line2 = "variableStep chrom="  # unknown
        with open(infile) as file:
            for line in file:
                if line.startswith("variableStep"):
                    line2 = line.rstrip()
                    break

        if self.combined_wig == True:
            (sites, data, files) = tnseq_tools.read_combined_wig(self.ctrldata[0])
        else:
            (data, sites) = tnseq_tools.CombinedWig.gather_wig_data(self.ctrldata)
        (data, factors) = norm_tools.normalize_data(data, self.normalization)

        print("writing", outputPath)
        file = open(outputPath, "w")
        file.write("# %s normalization of %s\n" % (self.normalization, infile))
        if self.combined_wig == True:
            for f in files:
                file.write("#File: %s\n" % f)
            for i in range(len(sites)):
                file.write(
                    "\t".join(
                        [str(sites[i])] + ["%0.1f" % x for x in list(data[..., i])]
                    )
                    + "\n"
                )
        else:
            file.write(line2 + "\n")
            for j in range(len(sites)):
                file.write("%s %s\n" % (sites[j], int(data[0, j])))
        file.close()

        self.finish()
        logging.log("Finished Normalization")

if __name__ == "__main__":

    (args, kwargs) = transit_tools.clean_args(sys.argv[1:])

    G = Norm.from_args(sys.argv[1:])
    G.Run()
