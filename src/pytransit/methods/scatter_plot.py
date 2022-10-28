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
from pytransit.generic_tools import csv, misc
from pytransit.specific_tools.transit_tools import wx, basename, HAS_R, FloatVector, DataFrame, StrVector
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components import samples_area, file_display, results_area, parameter_panel
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    name = "Scatter Plot"
    menu_name = name
    identifier = name.replace(" ",'')
    cli_name = identifier.lower()
    
    valid_cli_flags = [
        "-log",
        "-genes",
        "-cond" 
    ]

    usage_string = f"""
        usage:
            {console_tools.subcommand_prefix} {cli_name} <combined_wig> <metadata_file> <sample_id1> <sample_id2> <annotation_file> <output.png> [-genes -log]
            {console_tools.subcommand_prefix} {cli_name} <combined_wig> <metadata_file> -cond <condition1> <condition2> <annotation_file> <output.png> [-genes -log]
    """.replace("\n        ","\n     ")
    
    # usage: scatterplot <combined_wig> <metadata_file> <sample_id_or_condition1> <sample_id_or_condition2> <annotation_file> <output.png>
    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=6)
        
        combined_wig_path      = args[0]
        metadata_path          = args[1]
        sample_or_condition1   = args[2]
        sample_or_condition2   = args[3]
        annotation_path        = args[4]
        output_path            = args[5] # png file
        avg_by_conditions = "cond" in kwargs
        
        combined_wig = tnseq_tools.CombinedWig(
            main_path=combined_wig_path,
            metadata_path=metadata_path,
        )
        # 
        # filter either by condition or by sample
        # 
        if avg_by_conditions:
            combined_wig = combined_wig.with_only(condition_names=[sample_or_condition1, sample_or_condition2])
        else:
            combined_wig = combined_wig.with_only(wig_ids=[sample_or_condition1, sample_or_condition2])
        
        # save the data
        Method.output(
            combined_wig=combined_wig,
            annotation_path=annotation_path,
            output_path=output_path,
            normalization=kwargs["n"],
            gene_means="genes" in kwargs, # bool
            log_scale="log" in kwargs # bool
        )
        
    @gui.add_menu("Pre-Processing", cli_name)
    def on_menu_click(event):
        from pytransit.components import panel_helpers
        Method.value_getters = LazyDict()
        Method.by_condition = False
        with panel_helpers.NewPanel() as (panel, main_sizer):
            sample_ids = [x.id for x in gui.samples]
            Method.value_getters.combined_wig  = panel_helpers.combined_wig_filtered_by_sample_input(panel, main_sizer)
            Method.value_getters.gene_means    = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="average counts at the gene level", default_value=False, tooltip_text="if false, this shows the scatterplot of insertion counts at individual TA sites", widget_size=None)
            Method.value_getters.log_scale     = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="show axes on log scale", default_value=False, tooltip_text="show axes on log scale", widget_size=None)
            Method.value_getters.normalization = panel_helpers.create_normalization_input(panel, main_sizer)

            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=Method.from_gui)
    
    @staticmethod
    def from_gui(frame):
        arguments = LazyDict()
        
        # 
        # call all GUI getters, puts results into respective arguments
        # 
        for each_key, each_getter in Method.value_getters.items():
            try:
                arguments[each_key] = each_getter()
            except Exception as error:
                logging.error(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')

        # if plotting samples, (for now) only allow two samples
        # TODO: in future potentially allow a grid of scatterplots when more than two samples are selected
        if not Method.by_condition:
            assert len(arguments.combined_wig.samples) == 2, "Please select only two samples on the left"
        else:
            assert len(arguments.combined_wig.condition_names) == 2, "Please select only two conditions on the left"
        arguments.avg_by_conditions = Method.by_condition
        
        # 
        # ask for output path(s)
        # 
        arguments.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}.png",
            output_extensions='PNG file (*.png)|*.png;|\nAll files (*.*)|*.*',
        )
        arguments.normalization = "TTR" #TRI I should add a dropdown for this, but hard-code it for now
        # if user didn't select an output path
        if not arguments.output_path:
            return None

        Method.output(**arguments)

    @staticmethod
    def output(*, output_path, combined_wig_path=None, metadata_path=None, annotation_path=None, combined_wig=None, normalization=None, avg_by_conditions=None, gene_means=None, log_scale=None, n_terminus=None, c_terminus=None, disable_logging=False):
        # Defaults (even if argument directly provided as None)
        normalization     = normalization     if normalization     is not None else "TTR"
        avg_by_conditions = avg_by_conditions if avg_by_conditions is not None else False
        gene_means        = gene_means        if gene_means        is not None else False
        log_scale         = log_scale         if log_scale         is not None else False
        n_terminus        = n_terminus        if n_terminus        is not None else 0.0
        c_terminus        = c_terminus        if c_terminus        is not None else 0.0
        
        if combined_wig == None:
            combined_wig = tnseq_tools.CombinedWig(main_path=combined_wig_path, metadata_path=metadata_path, annotation_path=annotation_path)
        
        # TODO: in future potentially allow a grid of scatterplots when more than two samples are selected
        if not avg_by_conditions: assert len(combined_wig.samples) == 2, "Please use combined_wig.with_only(wig_ids=[ID1, ID2]) before calling scatter plot"
        else: assert len(combined_wig.condition_names) == 2, "Please use combined_wig.with_only(condition_names=[NAME1, NAME2]) before calling scatter plot"
        
        with transit_tools.TimerAndOutputs(method_name=Method.identifier, output_paths=[output_path], disable=disable_logging) as timer:
            # 
            # by gene or site
            # 
            if gene_means: 
                means, genes, labels = calc_gene_means(combined_wig, normalization, avg_by_conditions=avg_by_conditions, n_terminus=n_terminus, c_terminus=c_terminus)
                counts = means
            else:
                # average conditions if needed
                if avg_by_conditions:
                    combined_wig = combined_wig.averaged(by_conditions=True)
                counts = combined_wig.read_counts_array
                labels = combined_wig.wig_ids # wig_ids will be equivlent to condition names when averaged by conditions
            
            # 
            # plot
            # 
            sample1_name, sample2_name = labels
            sample_1_counts = counts[:,0]
            sample_2_counts = counts[:,1]
            import matplotlib.pyplot as plt
            if log_scale:
                plt.scatter(numpy.log10(sample_1_counts),numpy.log10(sample_2_counts))
                plt.xlabel("log10(%s)" % sample1_name)
                plt.ylabel("log10(%s)" % sample2_name)
            else:
                plt.scatter(sample_1_counts,sample_2_counts)
                plt.xlabel("%s" % sample1_name)
                plt.ylabel("%s" % sample2_name)
            if gene_means: plt.title("scatter plot of mean insertion counts for each gene")
            else: plt.title("scatter plot of insertion counts at individual TA sites")
            plt.savefig(output_path)
            plt.clf()