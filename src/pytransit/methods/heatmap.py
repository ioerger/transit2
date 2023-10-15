import collections
import datetime
import heapq
import math
import ntpath
import os
import random
import sys
import time

import numpy

from pytransit.components.spreadsheet import SpreadSheet
from pytransit.generic_tools import csv, informative_iterator, misc
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.globals import logging, cli, debugging_enabled, gui, root_folder
from pytransit.specific_tools import console_tools, gui_tools, norm_tools, tnseq_tools, transit_tools

@misc.singleton
class Method:
    name = "Heatmap"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name} - Perform {name} analysis"
    description = f"""Perform {name} analysis"""
    
    prev_menu_choice = None
    defaults = LazyDict(
        top_k=-1,
        q_value_threshold=0.05,
        low_mean_filter=5,
    )
    valid_cli_flags = [
        "-anova",
        "-zinb",
        "--qval",
        "--topk",
        "--low-mean-filter",
    ]
    usage_string = f"""
        Usage 1:
            {console_tools.subcommand_prefix} heatmap <anova_output> <heatmap.png> -anova [Optional Arguments]
        Usage 2:
            {console_tools.subcommand_prefix} heatmap <zinb_output> <heatmap.png> -zinb [Optional Arguments]
        
        Optional Arguments:
            --topk <int>            := number of results
            --qval <float>          := adjusted p value threshold. Default --qval 0.05
            --low-mean-filter <int> := Filter out genes with grand mean count (across all conditions) below this threshold
                                    (even if adjusted p-value < 0.05)
                                    Default --low-mean-filter 5
    """
    
    @gui.add_menu("Post-Processing", "ANOVA", "Heatmap")
    def on_menu_click(event):
        Method.prev_menu_choice = "anova"
        Method.define_panel(event)
    
    @gui.add_menu("Post-Processing", "ZINB", "Heatmap")
    def on_menu_click(event):
        Method.prev_menu_choice = "zinb"
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers, parameter_panel
        with panel_helpers.NewPanel() as (panel, main_sizer):
            parameter_panel.set_instructions(
                title_text=self.name,
                sub_text="",
                method_specific_instructions="""
                    The output of ANOVA or ZINB can be used to generate a heatmap that simultaneously clusters the significant genes and clusters the conditions, which is especially useful for shedding light on the relationships among the conditions apparent in the data.
                """.replace("\n                    ","\n"),
            )
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            
            self.value_getters = LazyDict()
            self.value_getters.input_path        = panel_helpers.create_file_input(  panel, main_sizer, button_label=f"Select {Method.prev_menu_choice} file", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*')
            self.value_getters.q_value_threshold = panel_helpers.create_float_getter(panel, main_sizer, label_text="Adj P Value Cutoff", default_value=Method.defaults.q_value_threshold, tooltip_text="Change adjusted p-value threshold for selecting genes")
            self.value_getters.top_k             = panel_helpers.create_int_getter(  panel, main_sizer, label_text="Top K",              default_value=Method.defaults.top_k,             tooltip_text="(-1 means all) Sometimes there are so many genes it is hard to see the heatmap top genes. This allows limiting to the top K genes (ranked by significance; adjusted p-value)")
            self.value_getters.low_mean_filter   = panel_helpers.create_float_getter(panel, main_sizer, label_text="Low Mean Filter",    default_value=Method.defaults.low_mean_filter,   tooltip_text="Filter out genes with grand mean count (across all conditions) below this threshold (even if adjusted p-value < 0.05)")
            
    @staticmethod
    def from_gui(frame):
        arguments = LazyDict()
        arguments.filetype = Method.prev_menu_choice # zinb or anova
        
        # 
        # call all GUI getters, puts results into respective arguments key-value
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
            default_file_name=f"{Method.cli_name}_output.png",
            output_extensions=transit_tools.result_output_extensions,
        )
        # if user didn't select an output path
        if not arguments.output_path:
            return None
        
        Method.load_from(**arguments)

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=2)
        
        if not kwargs["anova"] and not kwargs["zinb"]:
            logging.error(f"requires -anova or -zinb argument, see usage string below.\n{Method.usage_string}")
        
        Method.load_from(
            filetype="anova" if kwargs["anova"] else "zinb",
            input_path        = args[0],
            output_path       = args[1],
            q_value_threshold = kwargs["qval"],
            top_k             = kwargs["topk"],
            low_mean_filter   = kwargs["low-mean-filter"], # filter out genes with grandmean<5 by default
        )
        
    @staticmethod
    def load_from(filetype, input_path, output_path, top_k=None, q_value_threshold=None, low_mean_filter=None):
        if filetype == "anova":
            from pytransit.methods.anova import File as AnovaFile
            AnovaFile(path=input_path).create_heatmap(output_path, topk=top_k, qval=q_value_threshold, low_mean_filter=low_mean_filter)
        elif filetype == "zinb":
            from pytransit.methods.zinb import File as ZinbFile
            ZinbFile(path=input_path).create_heatmap(output_path, topk=top_k, qval=q_value_threshold, low_mean_filter=low_mean_filter)
        
    @staticmethod
    def output(column_names, formatted_rows, output_path, top_k=None, q_value_threshold=None, low_mean_filter=None):
        top_k              = int(top_k)               if top_k              != None  else Method.defaults.top_k
        q_value_threshold  = float(q_value_threshold) if q_value_threshold  != None  else Method.defaults.q_value_threshold
        low_mean_filter    = int(low_mean_filter)     if low_mean_filter    != None  else Method.defaults.low_mean_filter
        
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd
        
        # 
        # sort
        # 
        sorted_rows = sorted(formatted_rows, key=lambda row: row["q_value"])
        # 
        # filter by top_k
        # 
        slice_end = len(sorted_rows) if top_k == -1 else top_k
        sorted_rows = sorted_rows[:slice_end]
        # 
        # filter by significant
        # 
        significant_rows = [ each for each in sorted_rows if each["q_value"] < q_value_threshold ]
        # translation: olways have AT LEAST top_k elements in the sorted_rows
        if len(significant_rows) >= top_k:
            sorted_rows = significant_rows
        
        # 
        # apply low_mean_filter
        # 
        gene_names, lfc_s = [], []
        for each_row in sorted_rows:
            mean_of_means = round(numpy.mean(each_row["means"]), 1)
            if mean_of_means < low_mean_filter:
                print(f"""excluding {each_row["gene_name"]}, mean(means)={mean_of_means}""")
            else:
                gene_names.append(each_row["gene_name"])
                lfc_s.append(each_row['lfcs'])
        
        column_to_lfcs = pd.DataFrame({
            column_name : [ each_lfc[column_index] for each_lfc in lfc_s ]
                for column_index, column_name in enumerate(column_names)
        })

        basic_heatmap(column_to_lfcs, row_names=gene_names, output_path=output_path)
        
magic_number_300 = 300
magic_number_30 = 30
magic_number_15 = 15
def basic_heatmap(df, row_names, output_path):
    import matplotlib.pyplot as plt
    import seaborn as sns
    df.index = row_names
    
    number_of_columns = len(df.columns)
    number_of_rows = len(df.index)
    base_width  = magic_number_300+number_of_columns*magic_number_30
    base_height = magic_number_300+number_of_rows*magic_number_15
    scale = 1/plt.rcParams['figure.dpi'] 
    if os.path.isfile(output_path):
        os.remove(output_path) 
    plt.figure()
   
    clustermap_plot = sns.clustermap(
        df,
        figsize=(base_width*scale, base_height*scale),
        cmap=sns.color_palette('blend:Red,White,Blue', as_cmap=True),
        linewidths=.5,
        method="complete",
        metric="euclidean",
        center=0,
        yticklabels=True
    )
    x0, y0, cbar_width, cbar_height = clustermap_plot.cbar_pos
    clustermap_plot.ax_cbar.set_position([x0, 0.9, cbar_width/2, cbar_height])
    if os.path.exists(output_path):
        os.remove(output_path)
    plt.savefig(output_path, bbox_inches='tight')
