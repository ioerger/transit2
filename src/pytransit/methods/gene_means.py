import sys
import os
import time
import ntpath
import math
import random
import datetime
import collections
import heapq

import numpy

from pytransit.generic_tools import csv, misc, informative_iterator
from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.specific_tools.transit_tools import wx, basename
from pytransit.components.spreadsheet import SpreadSheet


@misc.singleton
class Method:
    name = "Gene Means"
    identifier  = name.replace(" ", "")
    cli_name    = name.replace(" ", "_").lower()
    menu_name   = f"Mean counts at the gene level ({name})"
    description = f"""Calculate mean counts at gene level."""
    
    valid_cli_flags = [
        "--n",  # normalization
        "--iN", # n_terminus
        "--iC", # c_terminus
        "-cond" # averages counts over replicates of each condition
    ]
    
    usage_string = f"""
        Usage:
            {console_tools.subcommand_prefix} {cli_name} <combined_wig> <metadata_file> <annotation_file> <output_file> [Optional Arguments]
        Optional Arguments:
            --n <string> :=  Normalization method. Default: --n TTR
            --iN <N>     :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: --iN 0
            --iC <N>     :=  Ignore TAs within given percentage (e.g. 5) of C terminus. Default: --iC 0
            -cond      :=  Averages counts over replicates of each condition
    """.replace("\n        ", "\n")
    
    @staticmethod
    @cli.add_command(cli_name)
    @cli.add_command("export", "gene_means")
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=4)
        
        # save the flags
        Method.output(
            combined_wig=tnseq_tools.CombinedWig.load(
                main_path=args[0],
                metadata_path=args[1],
                annotation_path=args[2],
            ),
            output_path=args[3],
            normalization=kwargs["n"],
            n_terminus=kwargs["iN"],
            c_terminus=kwargs["iC"],
            avg_by_conditions = "cond" in kwargs # boolean
        )
    
    @gui.add_menu("Pre-Processing", "Generate", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            self.value_getters = LazyDict()
            
            self.value_getters.n_terminus             = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus             = panel_helpers.create_c_terminus_input(panel, main_sizer)
            self.value_getters.normalization          = panel_helpers.create_normalization_input(panel, main_sizer)
            self.value_getters.avg_by_conditions      = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="Average counts over\nreplicates of each condition", default_value=False, tooltip_text="the output can contain gene means for each individual sample, or averaged by condition", widget_size=None)
            
            
    @staticmethod
    def from_gui(frame):
        arguments = LazyDict()
        
        # 
        # get global data
        # 
        arguments.combined_wig = gui.combined_wigs[-1]
        
        # 
        # call all GUI getters, puts results into respective key in arguments
        # 
        for each_key, each_getter in Method.value_getters.items():
            try:
                arguments[each_key] = each_getter()
            except Exception as error:
                logging.error(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
        
        # 
        # ask for output path(s)
        # 
        arguments.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}_output.tsv",
            output_extensions=transit_tools.result_output_extensions,
        )
        # if user didn't select an output path
        if not arguments.output_path:
            return None

        Method.output(**arguments)
    
    @staticmethod
    def calculate(combined_wig, avg_by_conditions=False, normalization="TTR", n_terminus=0, c_terminus=0):
        means, genes, labels = transit_tools.calc_gene_means(
            combined_wig=combined_wig,
            avg_by_conditions=avg_by_conditions,
            normalization=normalization,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
        )
        
        column_names = [
            "ORF",
            "Gene Name",
            "numTAs",
            *labels,
        ]

        rows = [ # expanded version of: for i in range(means.shape[0]): output.write("%s\n" % ('\t'.join([genes[i].orf,genes[i].name,genes[i].n]+["%0.1f" % x for x in means[i,:]])))
            [
                genes[row_index].orf,
                genes[row_index].name,
                genes[row_index].n,
                *[
                    "%0.1f" % x 
                        for x in means[row_index,:]
                ]
            ]
                for row_index in range(means.shape[0])
        ]
        return (column_names, rows), (means, genes, labels)
    
    @staticmethod
    def output(*, combined_wig, normalization=None, avg_by_conditions=None, output_path=None, n_terminus=None, c_terminus=None, disable_logging=False):
        # Defaults (even if argument directly provided as None)
        normalization     = normalization     if normalization     is not None else "TTR"
        avg_by_conditions = avg_by_conditions if avg_by_conditions is not None else False
        output_path       = output_path       if output_path       is not None else None
        n_terminus        = n_terminus        if n_terminus        is not None else 0.0
        c_terminus        = c_terminus        if c_terminus        is not None else 0.0

        from pytransit.specific_tools import stat_tools
        with transit_tools.TimerAndOutputs(method_name=Method.identifier, output_paths=[output_path], disable=disable_logging) as timer:
            (column_names, rows), (means, genes, labels) = Method.calculate(
                combined_wig=combined_wig,
                normalization=normalization,
                avg_by_conditions=avg_by_conditions,
                n_terminus=n_terminus,
                c_terminus=c_terminus,
            )

            logging.log(f"Writing output file: {output_path}")
            transit_tools.write_result(
                path=output_path, # path=None means write to STDOUT
                file_kind=Method.identifier,
                column_names=column_names,
                rows=rows,
                extra_info=dict(
                    time=timer.duration_in_seconds,
                ),
            )
        
@transit_tools.ResultsFile
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier,
                    heading=self.comments_string or misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[],
                ).Show(),
            })
        )
        
        # 
        # get data
        # 
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

