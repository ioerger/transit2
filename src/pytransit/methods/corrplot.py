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
    name = "Corrplot"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name} - Make correlation plot"
    description = f"""Make correlation plot"""
    
    transposons = [ "himar1" ] # not sure if this is right -- Jeff
    
    inputs = LazyDict(
        combined_wig=None,
        annotation_path=None,
        output_path=None,
        avg_by_conditions=False,
        normalization="TTR", #TRI hard-coded for now
    )
    
    valid_cli_flags = [ #TRI - consider adding these flags?
    ]
    # could add a flag for Adj P Value cutoff (or top n most signif genes)

    #TRI - should drop anova and zinb inputs, and instead take combined_wig or gene_means file (from export)
    #usage_string = """usage: {console_tools.subcommand_prefix} corrplot <gene_means> <output.png> [-anova|-zinb]""""
    usage_string = f"""usage: {console_tools.subcommand_prefix} corrplot <combined_wig> <annotation_file> <output.png> [-avg_by_conditions <metadata_file>]"""
    
    # 
    # CLI method
    # 
    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=3)
        
        # save the data
        Method.inputs.update(dict(
            combined_wig=tnseq_tools.CombinedWig(
                main_path=args[0],
                metadata_path=kwargs.get("avg_by_conditions",None),
            ),
            annotation_path=args[1],
            output_path=args[2], # png file,
            avg_by_conditions="avg_by_conditions" in kwargs, # bool
        ))
        
        Method.Run()
    
    # 
    # Button method
    # 
    @gui.add_wig_area_dropdown_option(name=name)
    def on_button_click(event):
        transit_tools.require_r_to_be_installed()
        
        selected_wig_fingerprints = [ each.fingerprint for each in gui.selected_samples or gui.samples ]
        Method.inputs.update(dict(
            combined_wig=gui.combined_wigs[-1].with_only(wig_fingerprints=selected_wig_fingerprints),
            annotation_path=gui.annotation_path,
            output_path=None,
            avg_by_conditions=False,
            normalization="TTR", #TRI hard-coded for now
        ))
        
        # 
        # validate
        # 
        transit_tools.validate_annotation(Method.inputs.annotation_path)
        assert len(selected_wig_fingerprints), "Please load at least one combined wig file"
        
        # 
        # Ask for inputs
        # 
        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"corrplot.png",
            output_extensions='PNG file (*.png)|*.png;|\nAll files (*.*)|*.*',
        )
        Method.Run()
    
    # 
    # Panel method
    # 
    @gui.add_menu("Pre-Processing", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        self.value_getters = LazyDict()
        with panel_helpers.NewPanel() as (panel, main_sizer):
            self.value_getters.avg_by_conditions = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="average counts by condition", default_value=False, tooltip_text="correlations among conditions (where counts are averaged among replicates of each condition) versus all individual samples", widget_size=None)
            self.value_getters.normalization     = panel_helpers.create_normalization_input(panel, main_sizer)

            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            
    @staticmethod
    def from_gui(frame):
        # 
        # get annotation
        # 
        Method.inputs.annotation_path = gui.annotation_path
        Method.inputs.combined_wig = gui.combined_wigs[-1] #TRI what if not defined? fail gracefully?
        
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
            default_file_name=f"corrplot.png",
            output_extensions='PNG file (*.png)|*.png;|\nAll files (*.*)|*.*',
        )
        
        # if user didn't select an output path
        if not Method.inputs.output_path:
            return None

        return Method
    
    # 
    # Run
    # 
    def Run(self):
        with gui_tools.nice_error_log:
            logging.log(f"Starting {Method.identifier} analysis")
            start_time = time.time()
            
            transit_tools.make_corrplot(
                combined_wig=self.inputs.combined_wig,
                normalization=self.inputs.normalization,
                annotation_path=self.inputs.annotation_path,
                avg_by_conditions=self.inputs.avg_by_conditions,
                output_path=self.inputs.output_path,
            )
            
            if gui.is_active:
                logging.log(f"Adding File: {self.inputs.output_path}")
                results_area.add(self.inputs.output_path)
            logging.log(f"Finished {Method.identifier} analysis in {time.time() - start_time:0.1f}sec")
        