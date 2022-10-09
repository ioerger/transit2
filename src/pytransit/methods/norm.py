from pytransit.components.parameter_panel import panel, progress_update
import sys
import os
import time

from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools, informative_iterator
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import misc
from pytransit.globals import gui, cli, root_folder, debugging_enabled

@misc.singleton
class Method:
    # FIXME: replace this Run() with the Run() from the full normalize.py
    name = "norm"
    description = "Method for normalizing datasets and outputting into CombinedWig file."
    usage_string = f"""
        {console_tools.subcommand_prefix} norm <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    
            Optional Arguments:
            -n <string>     :=  Normalization method. Default: -n TTR
    """.replace("\n        ","\n")

    
    inputs = LazyDict(
        ctrldata=None,
        normalization=None,
        annotation_path=None,
        output_file=None,
    )
    
    @cli.add_command("export", "norm")
    @staticmethod
    def from_args(args, kwargs):
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=3)

        Method.inputs.update(dict(
            ctrldata=args[0].split(","),
            annotation_path=args[1],
            output_path=args[2],
            normalization=kwargs.get("n", "TTR"),
        ))
        
        Method.Run()

    def Run(self):
        logging.log("Starting Normalization")
        
        transit_tools.convert_to_combined_wig(
            self.ctrldata,
            self.annotation_path,
            self.output_path,
            normchoice=self.normalization,
        )

        logging.log("Finished Normalization")
