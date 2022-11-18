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
from pytransit.specific_tools.transit_tools import wx, basename, r, globalenv, HAS_R, FloatVector, DataFrame, StrVector
from pytransit.components.spreadsheet import SpreadSheet


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
        transit_tools.require_r_to_be_installed()
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
        genenames = ["%s/%s" % (w[0], w[1]) for w in hits]
        hash = {}
        headers = [h.replace("Mean_", "") for h in headers]
        for i, col in enumerate(headers):
            hash[col] = FloatVector([x[i] for x in LFCs])
        df = DataFrame(hash)
        heatmapFunc = self.make_heatmap_r_func()
        heatmapFunc(df, StrVector(genenames), self.inputs.output_path)
        results_area.add(self.inputs.output_path)

    def make_heatmap_r_func(self):
        r("""
            make_heatmap = function(lfcs,genenames,outfilename) { 
                rownames(lfcs) = genenames
                suppressMessages(require(gplots))
                colors <- colorRampPalette(c("red", "white", "blue"))(n = 200)

                C = length(colnames(lfcs))
                R = length(rownames(lfcs))
                W = 300+C*30
                H = 300+R*15

                png(outfilename,width=W,height=H)
                #defaults are lwid=lhei=c(1.5,4)
                #heatmap.2(as.matrix(lfcs),col=colors,margin=c(12,12),lwid=c(2,6),lhei=c(0.1,2),trace="none",cexCol=1.4,cexRow=1.4,key=T) # make sure white=0
                #heatmap.2(as.matrix(lfcs),col=colors,margin=c(12,12),trace="none",cexCol=1.2,cexRow=1.2,key=T) # make sure white=0 # setting margins was causing failures, so remove it 8/22/21
                heatmap.2(as.matrix(lfcs),col=colors,margin=c(12,12),trace="none",cexCol=1.2,cexRow=1.2,key=T) # actually, margins was OK, so the problem must have been with lhei and lwid
                dev.off()
            }
        """)
        return globalenv["make_heatmap"]