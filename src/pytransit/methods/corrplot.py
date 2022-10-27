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
from pytransit.specific_tools.transit_tools import wx, r, basename, HAS_R, FloatVector, DataFrame, StrVector, globalenv
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    name = "Corrplot"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name} - Make correlation plot"
    description = f"""Make correlation plot"""
    
    transposons = [ "himar1" ] # not sure if this is right -- Jeff
    
    valid_cli_flags = [
        "-n",  # normalization
        "--avg_by_conditions",
        "-iN", # n_terminus
        "-iC", # c_terminus
        # could add a flag for Adj P Value cutoff (or top n most signif genes)
    ]

    # TODO: TRI - should drop anova and zinb defaults, and instead take combined_wig or gene_means file (from export)
    #usage_string = """usage: {console_tools.subcommand_prefix} corrplot <gene_means> <output.png> [-anova|-zinb]""""
    usage_string = f"""usage: {console_tools.subcommand_prefix} {cli_name} <combined_wig> <annotation_file> <output.png> [-avg_by_conditions <metadata_file>]"""
    
    # 
    # CLI method
    # 
    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=3)
        
        # map data to the core function
        Method.output(
            combined_wig=tnseq_tools.CombinedWig(
                main_path=args[0],
                metadata_path=kwargs.get("avg_by_conditions",None),
                annotation_path=args[1],
            ),
            normalization=kwargs["n"],
            n_terminus=kwargs["iN"],
            c_terminus=kwargs["iC"],
            avg_by_conditions="avg_by_conditions" in kwargs, # bool
            output_path=args[2],
            disable_logging=False,
        )
    
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
        arguments = LazyDict()
        
        # 
        # get global data
        # 
        arguments.combined_wig = gui.combined_wigs[-1] #TRI what if not defined? fail gracefully?
        
        # 
        # call all GUI getters, puts results into respective Method.defaults key-value
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
            default_file_name=f"corrplot.png",
            output_extensions='PNG file (*.png)|*.png;|\nAll files (*.*)|*.*',
        )
        
        # if user didn't select an output path
        if not arguments.output_path:
            return None
        
        # run the core function directly
        Method.output(**arguments)
    
    corrplot_r_function = None    
    @staticmethod
    def output(*, combined_wig, normalization=None, avg_by_conditions=None, output_path=None, n_terminus=None, c_terminus=None, disable_logging=False):
        # Defaults (even if argument directly provided as None)
        normalization     = normalization     if normalization     is not None else "TTR"
        avg_by_conditions = avg_by_conditions if avg_by_conditions is not None else False
        output_path       = output_path       if output_path       is not None else None
        n_terminus        = n_terminus        if n_terminus        is not None else 0.0
        c_terminus        = c_terminus        if c_terminus        is not None else 0.0
        
        from pytransit.methods.gene_means import Method as GeneMeansMethod
        with transit_tools.TimerAndOutputs(method_name=Method.identifier, output_paths=[output_path], disable=disable_logging,):
            import numpy
            # instantiate the corrplot_r_function if needed
            if not Method.corrplot_r_function:
                transit_tools.require_r_to_be_installed()
                r(""" # R function...
                    make_corrplot = function(means,headers,outfilename) { 
                        means = means[,headers] # put cols in correct order
                        suppressMessages(require(corrplot))
                        png(outfilename)
                        corrplot(cor(means))
                        dev.off()
                    }
                """)
                Method.corrplot_r_function = globalenv["make_corrplot"]
            
            _, (means, genes, labels) = GeneMeansMethod.calculate(combined_wig, normalization=normalization, avg_by_conditions=avg_by_conditions, n_terminus=n_terminus, c_terminus=c_terminus)
            print(f'''GeneMeansMethod.calculate: labels = {labels}''')
            
            position_hash = {}
            for i, col in enumerate(labels):
                position_hash[col] = FloatVector([x[i] for x in means])
            df = DataFrame(position_hash)  
            
            Method.corrplot_r_function(df, StrVector(labels), output_path ) # pass in headers to put cols in order, since df comes from dict