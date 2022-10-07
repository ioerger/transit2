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
from pytransit.tools.transit_tools import wx, pub, basename, FloatVector, DataFrame, StrVector, r, DataFrame, globalenv, IntVector, FloatVector, StrVector
from pytransit.universal_data import universal
from pytransit.components import file_display, results_area, parameter_panel, panel_helpers
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.panel_helpers import create_normalization_input, create_reference_condition_input, create_include_condition_list_input, create_exclude_condition_list_input, create_n_terminus_input, create_c_terminus_input, create_pseudocount_input, create_winsorize_input, create_alpha_input, create_button
from pytransit.interfaces import gui
command_name = sys.argv[0]

@misc.singleton
class Analysis:
    name = "Corrplot"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name} - Make correlation plot"
    description = f"""Make correlation plot"""
    
    transposons = [ "himar1" ] # not sure if this is right -- Jeff
    
    inputs = LazyDict(
        output_path=None,
        normalization="TTR",
        n_terminus=0.0,
        c_terminus=0.0,
    )
    
    valid_cli_flags = [ #TRI - consider adding these flags?
    ]
    # could add a flag for padj cutoff (or top n most signif genes)

    #TRI - should drop anova and zinb inputs, and instead take combined_wig or gene_means file (from export)
    #usage_string = """usage: python3 transit.py corrplot <gene_means> <output.png> [-anova|-zinb]""""
    usage_string = """usage: python3 transit.py corrplot <combined_wig> <annotation_file> <output.png> [-avg_by_conditions <metadata_file>]"""
    
    @gui.add_menu("Analysis", menu_name)
    def on_menu_click(event):
        Analysis.define_panel(event)
    
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
        Analysis.inputs.annotation_path =universal.annotation_path
        transit_tools.validate_annotation(Analysis.inputs.annotation_path)
        Analysis.inputs.combined_wig =universal.combined_wigs[0].main_path #TRI what if not defined? fail gracefully?
        Analysis.inputs.metadata =universal.combined_wigs[0].metadata.path
        
        # 
        # call all GUI getters, puts results into respective Analysis.inputs key-value
        # 
        for each_key, each_getter in Analysis.value_getters.items():
            try:
                Analysis.inputs[each_key] = each_getter()
            except Exception as error:
                logging.error(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
        #logging.log("included_conditions", Analysis.inputs.included_conditions)
        
        # 
        # ask for output path(s)
        # 
        Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"corrplot.png",
            output_extensions='PNG file (*.png)|*.png;|\nAll files (*.*)|*.*',
        )
        Analysis.inputs.normalization = "TTR" #TRI I should add a dropdown for this, but hard-code it for now
        # if user didn't select an output path
        if not Analysis.inputs.output_path:
            return None

        return Analysis

    @staticmethod
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Analysis.usage_string)
        console_tools.handle_unrecognized_flags(Analysis.valid_cli_flags, kwargs, Analysis.usage_string)
        console_tools.enforce_number_of_args(args, Analysis.usage_string, exactly=3)

        combined_wig    = args[0]
        annotation_path = args[1]
        output_path     = args[2] # png file
        metadata        = kwargs.get("avg_by_conditions",None) # optional
        avg_by_conditions = (metadata!=None) # boolean

        # save the data
        Analysis.inputs.update(dict(
            combined_wig = combined_wig,
            annotation_path = annotation_path,
            output_path = output_path,
            metadata = metadata,
            avg_by_conditions = avg_by_conditions, # bool
            normalization = "TTR", #TRI hard-coded for now
            #normalization=kwargs.get("n", Analysis.inputs.normalization), #TRI these might be useful for corrplot...
            #n_terminus=float(kwargs.get("iN", Analysis.inputs.n_terminus)),
            #c_terminus=float(kwargs.get("iC", Analysis.inputs.c_terminus)),
        ))
        
        return Analysis
        
    ##################################################

    def Run(self):
        with gui_tools.nice_error_log:
            logging.log(f"Starting {Analysis.identifier} analysis")
            start_time = time.time()
            
            transit_tools.make_corrplot(
                combined_wig=self.inputs.combined_wig,
                normalization=self.inputs.normalization,
                annotation_path=self.inputs.annotation_path,
                avg_by_conditions=self.inputs.avg_by_conditions,
                metadata=self.inputs.metadata,
                output_path=self.inputs.output_path,
            )
            
            if universal.interface=="gui":
                logging.log(f"Adding File: {self.inputs.output_path}")
                results_area.add(self.inputs.output_path)
            logging.log(f"Finished {Analysis.identifier} analysis in {time.time() - start_time:0.1f}sec")


Method = GUI = Analysis # for compatibility with older code/methods