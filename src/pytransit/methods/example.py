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

from pytransit.tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.basics.lazy_dict import LazyDict
import pytransit.basics.csv as csv
import pytransit.basics.misc as misc
from pytransit.tools.transit_tools import wx, pub, basename, HAS_R, FloatVector, DataFrame, StrVector, EOL
from pytransit.universal_data import universal
from pytransit.components import file_display, results_area, parameter_panel
from pytransit.components.spreadsheet import SpreadSheet
command_name = sys.argv[0]

name = "Example" # HANDLE_THIS

@misc.singleton
class Analysis:
    identifier  = name
    short_name  = name.lower()
    long_name   = name.upper()
    short_desc  = f"Perform {name} analysis"
    long_desc   = f"""Perform {name} analysis"""
    transposons = [ "himar1", "tn5" ]
    
    inputs = LazyDict(
        output_path=None,
        normalization="TTR",
        n_terminus=0.0,
        c_terminus=0.0,
        # HANDLE_THIS
    )
    
    valid_cli_flags = [
        "-n",  # normalization
        "-iN", # n_terminus
        "-iC", # c_terminus
        # HANDLE_THIS
    ]
    usage_string = f"""
        # HANDLE_THIS
        Usage: python3 transit.py {short_name} [Optional Arguments]
        Optional Arguments:
            -n <string>         :=  Normalization method. Default: -n TTR
            -iN <N> :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: -iN 0
            -iC <N> :=  Ignore TAs within given percentage (e.g. 5) of C terminus. Default: -iC 0
    """.replace("\n        ", "\n")
    
    
    wxobj = None
    panel = None
    
    def __init__(self, *args, **kwargs):
        self.full_name        = f"[{self.short_name}]  -  {self.short_desc}"
    
    def __str__(self):
        return f"""
            Analysis Method:
                Short Name:  {self.short_name}
                Long Name:   {self.long_name}
                Short Desc:  {self.short_desc}
                Long Desc:   {self.long_desc}
        """.replace('\n            ','\n').strip()
    
    def __repr__(self): return f"{self.inputs}"
    def __call__(self): return self
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (self.panel, main_sizer):
            self.value_getters = LazyDict()
            # panel_helpers.create_float_getter(self.panel, main_sizer, label_text="", default_value=0, tooltip_text="")
            # panel_helpers.create_int_getter(self.panel, main_sizer, label_text="", default_value=0, tooltip_text="")
            # panel_helpers.create_file_input(self.panel, main_sizer, button_label="", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*')
            # panel_helpers.create_choice_input(self.panel, main_sizer, label="", options=[], default_option=None, tooltip_text="")
            # panel_helpers.create_text_box_getter(self.panel, main_sizer, label_text="", default_value="", tooltip_text="", label_size=None, widget_size=None,)
            # panel_helpers.create_check_box_getter(self.panel, main_sizer, label_text="", default_value=False, tooltip_text="", widget_size=None)
            # @panel_helpers.create_button(self.panel, main_sizer, label="")
            # def when_button_clicked(event):
            #     print("do stuff")
            
            self.value_getters.n_terminus             = panel_helpers.create_n_terminus_input(self.panel, main_sizer)
            self.value_getters.c_terminus             = panel_helpers.create_c_terminus_input(self.panel, main_sizer)
            self.value_getters.normalization          = panel_helpers.create_normalization_input(self.panel, main_sizer)
            
            panel_helpers.create_run_button(self.panel, main_sizer, from_gui_function=self.from_gui)
            

    @staticmethod
    def from_gui(frame):
        # 
        # global data
        # 
        # HANDLE_THIS
        universal.interface # "gui" or "console"
        universal.frame # self.wxobj equivalent
        universal.busy_running_method # Boolean, is true when any .Run() is started but not finished
        universal # I would like to flatten this (remove the .session_data namespace) but its low priority
        universal.annotation_path # string, may need to become a list of strings
        universal.conditions # list of Condition objects
        universal.conditions[0].name # string
        universal.conditions[0].extra_data # dict (currently unused, but would show up as columns in the condition GUI table)
        universal.combined_wigs # list of CombinedWig objects
        universal.combined_wigs[0].main_path
        universal.combined_wigs[0].metadata_path # to get all these it would be [ each.metadata_path for each in universal.combined_wigs ]
        universal.combined_wigs[0].samples # list of Wig objects
        universal.combined_wigs[0].samples[0].id # id from the metadata file
        universal.combined_wigs[0].samples[0].fingerprint # the "File" column from the metadata 
        universal.combined_wigs[0].samples[0].condition_names # a list of strings
        universal.combined_wigs[0].samples[0].positions # list of ints
        universal.combined_wigs[0].samples[0].insertion_counts # list of numbers
        universal.combined_wigs[0].samples[0].rows # each element is always [position_number, insertion_count]
        universal.combined_wigs[0].samples[0].column_index # int (column inside combined wig)
        universal.combined_wigs[0].samples[0].extra_data.count
        universal.combined_wigs[0].samples[0].extra_data.sum
        universal.combined_wigs[0].samples[0].extra_data.non_zero_mean
        universal.combined_wigs[0].samples[0].extra_data.non_zero_median
        universal.combined_wigs[0].samples[0].extra_data.density
        universal.combined_wigs[0].samples[0].extra_data.mean
        universal.combined_wigs[0].samples[0].extra_data.max
        universal.combined_wigs[0].samples[0].extra_data.skew
        universal.combined_wigs[0].samples[0].extra_data.kurtosis
        universal.combined_wigs[0].metadata # CombinedWigMetadata object
        universal.combined_wigs[0].metadata.path
        universal.combined_wigs[0].metadata.headers
        universal.combined_wigs[0].metadata.rows
        universal.combined_wigs[0].metadata.conditions
        universal.combined_wigs[0].metadata.condition_for(wig_fingerprint) # will need to change to "conditions" instead of "condition"
        universal.combined_wigs[0].metadata.condition_for(wig_id) # will need to change to "conditions" instead of "condition"
        universal.combined_wigs[0].metadata.id_for(wig_fingerprint)
        universal.combined_wigs[0].metadata.fingerprints_for(condition_name)
        universal.combined_wigs[0].rows # equivalent to the CSV rows of .comwig file; a list of lists, can contain numbers and strings
        
        # 
        # get annotation
        # 
        # HANDLE_THIS
        Analysis.inputs.annotation_path = universal.annotation_path
        transit_tools.validate_annotation(Analysis.inputs.annotation_path)
        
        # 
        # call all GUI getters, puts results into respective Analysis.inputs key-value
        # 
        for each_key, each_getter in Analysis.value_getters.items():
            try:
                Analysis.inputs[each_key] = each_getter()
            except Exception as error:
                logging.error(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
        logging.log("included_conditions", Analysis.inputs.included_conditions)
        
        # 
        # ask for output path(s)
        # 
        Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Analysis.short_name}_output.dat",
            output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
        )
        # if user didn't select an output path
        if not Analysis.inputs.output_path:
            return None

        return Analysis

    @staticmethod
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Analysis.usage_string)
        console_tools.handle_unrecognized_flags(Analysis.valid_cli_flags, kwargs, Analysis.usage_string)

        # save the data
        Analysis.inputs.update(dict(
            output_path=args[0],
            normalization=kwargs.get("n", Analysis.inputs.normalization),
            n_terminus=float(kwargs.get("iN", Analysis.inputs.n_terminus)),
            c_terminus=float(kwargs.get("iC", Analysis.inputs.c_terminus)),
            # HANDLE_THIS
        ))
        
        return Analysis
        
    def Run(self):
        from pytransit.tools import stat_tools
        logging.log(f"Starting {Analysis.identifier} analysis")
        start_time = time.time()
        
        # 
        # process data
        # 
        if True:
            rows, summary_info = stat_tools.{analysis_name}(**self.inputs) # HANDLE_THIS
        
        # 
        # write output
        # 
        if True:
            logging.log(f"Adding File: {self.inputs.output_path}")
            # 
            # write to file
            # 
            transit_tools.write_result(
                path=self.inputs.output_path, # path=None means write to STDOUT
                file_kind=Analysis.identifier,
                rows=rows,
                column_names=[
                    # HANDLE_THIS
                ],
                extra_info=dict(
                    stats=dict(summary_info), # HANDLE_THIS
                    parameters=self.inputs,
                ),
            )
            logging.log(f"Finished {Analysis.identifier} analysis in {time.time() - start_time:0.1f}sec")
        results_area.add(self.inputs.output_path)

@transit_tools.ResultsFile
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Analysis.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Analysis.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Analysis.identifier,
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[
                        # HANDLE_THIS
                    ],
                ).Show(),
            })
        )
        
        # 
        # get column names
        # 
        self.column_names, self.rows, self.extra_data = tnseq_tools.read_results_file(self.path)
        self.values_for_result_table.update(self.extra_data.get("parameters", {}))
        
        # 
        # get summary stats
        #
        self.values_for_result_table.update({
            # HANDLE_THIS (additional summary_info for results table)
            # examples:
                # f"Gene Count": len(self.rows),
                # f"Padj<{Analysis.significance_threshold}": len([
                #     1 for each in self.rows
                #         if each.get("Padj", 0) < Analysis.significance_threshold 
                # ]),
        })
    
    def __str__(self):
        return f"""
            File for {Analysis.short_name}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

Analysis.filetypes = [ ResultFileType1, ] # HANDLE_THIS
Method = GUI = Analysis # for compatibility with older code/methods