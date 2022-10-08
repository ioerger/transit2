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
from pytransit.tools import transit_tools, tnseq_tools, norm_tools, stat_tools, console_tools
from pytransit.globals import gui, cli, root_folder, debugging_enabled


############# Description ##################

short_name = "norm"
long_name = "Normalization"
short_desc = "Normalization method"
long_desc = "Method for normalizing datasets and outputting into CombinedWig file."
transposons = ["himar1", "tn5"]
columns = ["Position", "Reads", "Genes"]


############# Analysis Method ##############


class Export(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            transposons,
            NormMethod,
            NormGUI,
            [NormFile],
        )


################## FILE ###################


class NormFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#CombinedWig", columns)

    def get_header(self, path):
        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
        return text


################# GUI ##################


class NormGUI(base.AnalysisGUI):
    def __init__(self):
        base.AnalysisGUI.__init__(self)


########## METHOD #######################


class NormMethod(base.SingleConditionMethod):
    """   
    Norm
 
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

    @cli.add_command("export", "norm")
    @staticmethod
    def from_args(args, kwargs):
        from pytransit.tools.console_tools import InvalidArgumentException
        if len(args) < 3:
            raise InvalidArgumentException("Must provide all necessary arguments")

        ctrldata = args[0].split(",")
        annotation_path = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        replicates = "Sum"
        normalization = kwargs.get("n", "TTR")
        LOESS = False
        ignore_codon = True
        n_terminus = 0.0
        c_terminus = 0.0

        return NormMethod(
            ctrldata,
            annotation_path,
            output_file,
            replicates,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
        )

    def Run(self):

        logging.log("Starting Normalization")
        start_time = time.time()
        outputPath = self.output.name
        # Normalize Data
        transit_tools.convert_to_combined_wig(
            self.ctrldata,
            self.annotation_path,
            outputPath,
            normchoice=self.normalization,
        )

        self.finish()
        logging.log("Finished Normalization")

    usage_string = f"""
{console_tools.subcommand_prefix} norm <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    
        Optional Arguments:
        -n <string>     :=  Normalization method. Default: -n TTR
        """