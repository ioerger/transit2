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

from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
import pytransit.generic_tools.csv as csv
import pytransit.generic_tools.misc as misc
from pytransit.specific_tools.transit_tools import wx, pub, basename, FloatVector, DataFrame, StrVector, r, DataFrame, globalenv, IntVector, FloatVector, StrVector
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components import file_display, results_area, parameter_panel, panel_helpers
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.panel_helpers import create_normalization_input, create_reference_condition_input, create_include_condition_list_input, create_exclude_condition_list_input, create_n_terminus_input, create_c_terminus_input, create_pseudocount_input, create_winsorize_input, create_alpha_input, create_button

@misc.singleton
class Method:
    name = "Scatterplot"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name}"
    description = f"""Make scatter plot of insertion counts between samples, averaged at the gene level"""
    
    transposons = [ "himar1" ] # not sure if this is right -- Jeff
    
    inputs = LazyDict(
        output_path=None,
        normalization="TTR",
        n_terminus=0.0,
        c_terminus=0.0,
    )
    
    valid_cli_flags = [ #TRI - consider adding these flags?
    ]

    usage_string = f"""usage: {console_tools.subcommand_prefix} scatterplot <combined_wig> <metadata_file> <sample_id_or_condition1> <sample_id_or_condition2> <annotation_file> <output.png>"""
    
    @gui.add_menu("Method", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        self.value_getters = LazyDict()
        with panel_helpers.NewPanel() as (panel, main_sizer):
            self.value_getters.avg_by_conditions = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="average counts by condition", default_value=False, tooltip_text="correlations among conditions (where counts are averaged among replicates of each condition) versus all individual samples", widget_size=None)

            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            
    @staticmethod
    def from_gui(frame):
        # 
        # get annotation
        # 
        Method.inputs.annotation_path =gui.annotation_path
        transit_tools.validate_annotation(Method.inputs.annotation_path)
        Method.inputs.combined_wig =gui.combined_wigs[0].main_path #TRI what if not defined? fail gracefully?
        Method.inputs.metadata =gui.combined_wigs[0].metadata.path
        
        # 
        # call all GUI getters, puts results into respective Method.inputs key-value
        # 
        for each_key, each_getter in Method.value_getters.items():
            try:
                Method.inputs[each_key] = each_getter()
            except Exception as error:
                logging.error(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
        #logging.log("included_conditions", Method.inputs.included_conditions)
        
        # 
        # ask for output path(s)
        # 
        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"corrplot.png",
            output_extensions='PNG file (*.png)|*.png;|\nAll files (*.*)|*.*',
        )
        Method.inputs.normalization = "TTR" #TRI I should add a dropdown for this, but hard-code it for now
        # if user didn't select an output path
        if not Method.inputs.output_path:
            return None

        return Method


    # usage: scatterplot <combined_wig> <metadata_file> <sample_id_or_condition1> <sample_id_or_condition2> <annotation_file> <output.png>

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=6)

        combined_wig    = args[0]
        metadata        = args[1]
        sample1         = args[2]
        sample2         = args[3]
        annotation_path = args[4]
        output_path     = args[5] # png file

        # save the data
        Method.inputs.update(dict(
            combined_wig = combined_wig,
            metadata = metadata,
            annotation_path = annotation_path,
            sample1 = sample1,
            sample2 = sample2,
            output_path = output_path,
            normalization=kwargs.get("n", "TTR")
        ))
        
        Method.Run()
        
    ##################################################

    def Run(self):
      print("here I am")
      if False:
        with gui_tools.nice_error_log:
            logging.log(f"Starting {Method.identifier} analysis")
            start_time = time.time()
            
            self.make_scatterplot(
                combined_wig=self.inputs.combined_wig,
                normalization=self.inputs.normalization,
                annotation_path=self.inputs.annotation_path,
                avg_by_conditions=self.inputs.avg_by_conditions,
                metadata=self.inputs.metadata,
                output_path=self.inputs.output_path,
                sample1=self.inputs.sample1,
                sample2=self.inputs.sample2,
            )
            
            if gui.is_active:
                logging.log(f"Adding File: {self.inputs.output_path}")
                results_area.add(self.inputs.output_path)
            logging.log(f"Finished {Method.identifier} analysis in {time.time() - start_time:0.1f}sec")


