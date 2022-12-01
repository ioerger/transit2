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

from pytransit.components import file_display, parameter_panel, results_area, samples_area
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.generic_tools import csv, informative_iterator, misc
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.globals import cli, debugging_enabled, gui, root_folder
from pytransit.specific_tools import console_tools, gui_tools, logging, norm_tools, tnseq_tools, transit_tools


@misc.singleton
class Method:
    name = "Heatmap"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name} - Perform {name} analysis"
    description = f"""Perform {name} analysis"""
    
    inputs = LazyDict(
        filetype=None, # "anova" or "zinb"
        input_path=None,
        output_path=None,
        adj_p_value=0.05,
        top_k=-1,
        low_mean_filter=5,
    )
    valid_cli_flags = [
        "--anova",
        "--zinb",
        "-qval",
        "-topk",
        "-low-mean-filter",
    ]
    usage_string = f"""
        Usage 1:
            {console_tools.subcommand_prefix} heatmap <anova_output> <heatmap.png> --anova [Optional Arguments]
        Usage 2:
            {console_tools.subcommand_prefix} heatmap <zinb_output> <heatmap.png> --zinb [Optional Arguments]
        
        Optional Arguments:
            -topk <int>            := number of results
            -qval <float>          := adjusted p value threshold. Default -qval 0.05
            -low-mean-filter <int> := Filter out genes with grand mean count (across all conditions) below this threshold
                                    (even if adjusted p-value < 0.05)
                                    Default -low-mean-filter 5
    """
    
    @gui.add_menu("Post-Processing", "ANOVA", "Heatmap")
    def on_menu_click(event):
        Method.inputs.filetype = "anova"
        Method.define_panel(event)
    
    @gui.add_menu("Post-Processing", "ZINB", "Heatmap")
    def on_menu_click(event):
        Method.inputs.filetype = "zinb"
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            parameter_panel.set_instructions(
                title_text=self.name,
                sub_text="",
                method_specific_instructions="""
                    The output of ANOVA or ZINB can be used to generate a heatmap that simultaneously clusters the significant genes and clusters the conditions, which is especially useful for shedding light on the relationships among the conditions apparent in the data.
                    Note: The heatmap command calls R, which must be installed on your system, and relies on the 'gplots' R package.
                """.replace("\n                    ","\n"),
            )
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            
            self.value_getters = LazyDict()
            self.value_getters.input_path      = panel_helpers.create_file_input(  panel, main_sizer, button_label=f"Select {Method.inputs.filetype} file", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*')
            self.value_getters.adj_p_value     = panel_helpers.create_float_getter(panel, main_sizer, label_text="Adj P Value",     default_value=Method.inputs.adj_p_value,     tooltip_text="Change adjusted p-value threshold for selecting genes")
            self.value_getters.top_k           = panel_helpers.create_int_getter(  panel, main_sizer, label_text="Top K",           default_value=Method.inputs.top_k,           tooltip_text="Select top k genes ranked by significance (adjusted pval)")
            self.value_getters.low_mean_filter = panel_helpers.create_float_getter(panel, main_sizer, label_text="Low Mean Filter", default_value=Method.inputs.low_mean_filter, tooltip_text="Filter out genes with grand mean count (across all conditions) below this threshold (even if adjusted p-value < 0.05)")
            
            
    @staticmethod
    def from_gui(frame):
        # 
        # get annotation
        # 
        Method.inputs.annotation_path = gui.annotation_path
        
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
            default_file_name=f"{Method.cli_name}_output.png",
            output_extensions=transit_tools.result_output_extensions,
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
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=2)
        
        Method.inputs.filetype = None
        if kwargs.get("anova", False):
            Method.inputs.filetype = "anova"
        elif kwargs.get("zinb", False):
            Method.inputs.filetype = "zinb"
        else:
            logging.error(f"requires --anova or --zinb argument, see usage string below.\n{Method.usage_string}")
        
        Method.inputs.input_path = args[0]
        Method.inputs.output_path = args[1]
        Method.inputs.adj_p_value = float(kwargs.get("qval", 0.05))
        Method.inputs.top_k = int(kwargs.get("topk", -1))
        Method.inputs.low_mean_filter = int(
            kwargs.get("low-mean-filter", 5)
        )  # filter out genes with grandmean<5 by default
        
        Method.Run()
        
    def Run(self):
        if self.inputs.filetype != "anova" and self.inputs.filetype != "zinb":
            logging.error("filetype not recognized: %s" % self.inputs.filetype)

        headers = None
        data, hits = [], []
        n = -1  # number of conditions

        with open(self.inputs.input_path) as file:
            for line in file:
                w = line.rstrip().split("\t")
                if line[0] == "#" or ("P Value" in line):  # check for 'pval' for backwards compatibility
                    headers = w
                    continue  # keep last comment line as headers
                # assume first non-comment line is header
                if n == -1:
                    # ANOVA header line has names of conditions, organized as 3+2*n+3 (2 groups (means, LFCs) X n conditions)
                    # ZINB header line has names of conditions, organized as 3+4*n+3 (4 groups X n conditions)
                    if self.inputs.filetype == "anova":
                        n = int((len(w) - 6) / 2)
                    elif self.inputs.filetype == "zinb":
                        n = int((len(headers) - 6) / 4)
                    headers = headers[3 : 3 + n]
                    headers = [x.replace("Mean_", "") for x in headers]
                else:
                    means = [
                        float(x) for x in w[3 : 3 + n]
                    ]  # take just the columns of means
                    lfcs = [
                        float(x) for x in w[3 + n : 3 + n + n]
                    ]  # take just the columns of LFCs
                    qval = float(w[-2])
                    data.append((w, means, lfcs, qval))

        data.sort(key=lambda x: x[-1])
        hits, LFCs = [], []
        for k, (w, means, lfcs, qval) in enumerate(data):
            if (self.inputs.top_k == -1 and qval < self.inputs.adj_p_value) or (
                self.inputs.top_k != -1 and k < self.inputs.top_k
            ):
                mm = round(numpy.mean(means), 1)
                if mm < self.inputs.low_mean_filter:
                    print("excluding %s/%s, mean(means)=%s" % (w[0], w[1], mm))
                else:
                    hits.append(w)
                    LFCs.append(lfcs)

        print("heatmap based on %s genes" % len(hits))
        gene_names = ["%s/%s" % (w[0], w[1]) for w in hits]
        
        with transit_tools.TimerAndOutputs(method_name=Method.identifier, output_paths=[self.inputs.output_path],):
            import numpy as np
            import pandas as pd

            hash = {}
            headers = [h.replace("Mean_", "") for h in headers]
            for i, col in enumerate(headers):
                hash[col] = FloatVector([x[i] for x in LFCs])
            df = pd.DataFrame.from_dict(hash, orient="columns") 
           
            make_heatmap(df, gene_names, self.inputs.output_path)

    # @staticmethod
    # def output():
    #     # 
    #     # sort
    #     # 
    #     sorted_rows = sorted(formatted_rows, key=lambda row: row["q_value"])
    #     # 
    #     # filter by top_k
    #     # 
    #     slice_end = len(sorted_rows) if top_k == -1 else top_k
    #     sorted_rows = sorted_rows[:slice_end]
    #     # 
    #     # filter by significant
    #     # 
    #     significant_rows = [ each for each in sorted_rows if each["q_value"] < qval_threshold ]
    #     # translation: olways have AT LEAST top_k elements in the sorted_rows
    #     if len(significant_rows) >= top_k:
    #         sorted_rows = significant_rows
        
    #     # 
    #     # apply low_mean_filter
    #     # 
    #     gene_names, lfc_s = [], []
    #     for each_row in sorted_rows:
    #         mean_of_means = round(numpy.mean(each_row["means"]), 1)
    #         if mean_of_means < low_mean_filter:
    #             print(f"""excluding {each_row["gene_name"]}, mean(means)={mean_of_means}""")
    #         else:
    #             gene_names.append(each_row["gene_name"])
    #             lfc_s.append(each_row['lfcs'])

magic_number_300 = 300
magic_number_30 = 30
magic_number_15 = 15
def make_heatmap(df, gene_names, output_path):
    import matplotlib.pyplot as plt
    import seaborn as sns
    df.index = gene_names
    
    number_of_columns = len(df.columns)
    number_of_rows = len(df.index)
    base_width  = magic_number_300+number_of_columns*magic_number_30
    base_height = magic_number_300+number_of_rows*magic_number_15
    scale = 1/plt.rcParams['figure.dpi'] 
    if os.path.isfile(output_path):
        os.remove(output_path) 
    plt.figure()
   
    g = sns.clustermap(
        df,
        figsize=(base_width*scale, base_height*scale),
        cmap=sns.color_palette('blend:Red,White,Blue', as_cmap=True),
        linewidths=.5,
        method="complete",
        metric="euclidean",
        center=0,
    )
    x0, y0, cbar_width, cbar_height = g.cbar_pos
    g.ax_cbar.set_position([x0, 0.9, cbar_width/2, cbar_height])
    if os.path.exists(output_path):
        os.remove(output_path)
    plt.savefig(output_path, bbox_inches='tight')


def save_significance_heatmap(column_names, formatted_rows, output_path, top_k=-1, qval_threshold=0.05, low_mean_filter=5):
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
    significant_rows = [ each for each in sorted_rows if each["q_value"] < qval_threshold ]
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
    
    print(f"heatmap based on {len(gene_names)} genes")
    column_to_lfcs = pd.DataFrame({
        column_name : [ each_lfc[column_index] for each_lfc in lfc_s ]
            for column_index, column_name in enumerate(column_names)
    })
    
    column_to_lfcs.index = gene_names
    
    number_of_columns = len(column_to_lfcs.columns)
    number_of_rows = len(column_to_lfcs.index)
    base_width  = magic_number_300+number_of_columns*magic_number_30
    base_height = magic_number_300+number_of_rows*magic_number_15
    scale = 1/plt.rcParams['figure.dpi'] 
    if os.path.isfile(output_path):
        os.remove(output_path) 
    plt.figure()
   
    g = sns.clustermap(
        column_to_lfcs,
        figsize=(base_width*scale, base_height*scale),
        cmap=sns.color_palette('blend:Red,White,Blue', as_cmap=True),
        linewidths=.5,
        method="complete",
        metric="euclidean",
        center=0,
    )
    x0, y0, cbar_width, cbar_height = g.cbar_pos
    g.ax_cbar.set_position([x0, 0.9, cbar_width/2, cbar_height])
    if os.path.exists(output_path):
        os.remove(output_path)
    plt.savefig(output_path, bbox_inches='tight')