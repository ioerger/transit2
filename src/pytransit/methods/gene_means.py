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
from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components import samples_area, results_area, parameter_panel, file_display

from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.specific_tools.transit_tools import wx, basename, HAS_R, FloatVector, DataFrame, StrVector
from pytransit.components.spreadsheet import SpreadSheet


@misc.singleton
class Method:
    name = "Gene Means"
    identifier  = name.replace(" ", "")
    cli_name    = name.replace(" ", "_").lower()
    menu_name   = f"{name} - calculate mean counts at gene level"
    description = f"""Calculate mean counts at gene level."""
    
    inputs = LazyDict(
        combined_wig=None,
        metadata=None,
        annotation_path=None,
        normalization="TTR",
        condition_avg=None,
        output_path=None,
        
        n_terminus=0.0, # TODO: these don't seem to be used in the Run funciton --Jeff
        c_terminus=0.0, # TODO: these don't seem to be used in the Run funciton --Jeff
    )
    
    valid_cli_flags = [
        "-n",  # normalization
        "-iN", # n_terminus
        "-iC", # c_terminus
        "-cond" # averages counts over replicates of each condition
    ]
    
    usage_string = f"""
        Usage: {console_tools.subcommand_prefix} {cli_name} <combined_wig> <metadata> <prot_table> <output_file> [Optional Arguments]
        Optional Arguments:
            -n <string>         :=  Normalization method. Default: -n TTR
            -iN <N> :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: -iN 0
            -iC <N> :=  Ignore TAs within given percentage (e.g. 5) of C terminus. Default: -iC 0
            -cond   :=  Averages counts over replicates of each condition
    """.replace("\n        ", "\n")
    
    @gui.add_menu("Pre-Processing", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            self.value_getters = LazyDict()
            
            self.value_getters.n_terminus             = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus             = panel_helpers.create_c_terminus_input(panel, main_sizer)
            self.value_getters.normalization          = panel_helpers.create_normalization_input(panel, main_sizer)
            self.value_getters.condition_avg          = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="average counts over replicates of each condition", default_value=False, tooltip_text="the output can contain gene means for each individual sample, or averaged by condition", widget_size=None)
            
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            
    @staticmethod
    def from_gui(frame):
        # 
        # get annotation
        # 
        Method.inputs.combined_wig = gui.combined_wigs[-1] # what if user wants to change normalization or terminus trimming?
        Method.inputs.annotation_path = gui.annotation_path
        transit_tools.validate_annotation(gui.annotation_path)
        
        # 
        # call all GUI getters, puts results into respective Method.inputs key-value
        # 
        for each_key, each_getter in Method.value_getters.items():
            try:
                Method.inputs[each_key] = each_getter()
            except Exception as error:
                logging.error(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
        
        # 
        # ask for output path(s)
        # 
        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}_output.csv",
            output_extensions='Common output extensions (*.csv,*.txt,*.dat,*.out)|*.csv;*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
        )
        # if user didn't select an output path
        if not Method.inputs.output_path:
            return None

        return Method

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=4)
        
        # save the flags
        Method.inputs.update(dict(
            combined_wig_path=args[0],
            metadata_path=args[1],
            annotation_path=args[2],
            output_path=args[3],
            normalization=kwargs.get("n", Method.inputs.normalization),
            n_terminus=float(kwargs.get("iN", Method.inputs.n_terminus)),
            c_terminus=float(kwargs.get("iC", Method.inputs.c_terminus)),
            condition_avg = "cond" in kwargs, # boolean
            ) )
        Method.inputs.update(dict(
            combined_wig=tnseq_tools.CombinedWig(
                main_path=Method.inputs.combined_wig_path,
                metadata_path=Method.inputs.metadata_path,
                comments=None,
                extra_data=None
            ) ) ) 
        Method.Run()
        
    def Run(self):
        from pytransit.specific_tools import stat_tools
        logging.log(f"Starting {Method.identifier} analysis")
        start_time = time.time()
        
        # 
        # process data
        # 
        means, genes, labels = transit_tools.calc_gene_means(
            combined_wig_path=self.inputs.combined_wig_path,
            metadata_path=self.inputs.metadata_path,
            combined_wig = self.inputs.combined_wig,
            annotation_path=self.inputs.annotation_path,
            normalization=self.inputs.normalization,
            n_terminus=self.inputs.n_terminus,
            c_terminus=self.inputs.c_terminus,
            avg_by_conditions=self.inputs.condition_avg,
        ) 
        
        #
        # write output
        #
        rows = [ # expanded version of: for i in range(means.shape[0]): output.write("%s\n" % ('\t'.join([genes[i].orf,genes[i].name]+["%0.1f" % x for x in means[i,:]])))
            [
                genes[row_index].orf,
                genes[row_index].name,
                *[
                    "%0.1f" % x 
                        for x in means[row_index,:]
                ]
            ]
                for row_index in range(means.shape[0])
        ]
        logging.log(f"Writing output file: {self.inputs.output_path}")
        transit_tools.write_result(
            path=self.inputs.output_path, # path=None means write to STDOUT
            file_kind=Method.identifier,
            column_names=[
                "ORF",
                "Gene Name",
                *labels,
            ],
            rows=rows,
            extra_info=dict(
                time=(time.time() - start_time),
            ),
        )
        
        results_area.add(self.inputs.output_path)
        logging.log(f"Finished {Method.identifier} analysis in {time.time() - start_time:0.1f}sec")

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

