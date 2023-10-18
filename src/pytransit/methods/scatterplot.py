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

from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import csv, misc
from pytransit.specific_tools.transit_tools import calc_gene_means, wx, basename
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    name = "Scatter Plot"
    menu_name = f"{name}"
    identifier = name.replace(" ",'')
    cli_name = identifier.lower()
    
    valid_cli_flags = [
        "-log",
        "-genes",
        "--cond",
        "--samp",
    ]

    usage_string = f"""
        Usage:
            {console_tools.subcommand_prefix} {cli_name} <combined_wig_file> <metadata_file> <annotation_file> --samp <comma-separated sample ID's> <output.png> [Optional Arguments]
            Optional Arguments:
                -genes         := If set, this shows the scatterplot of mean insertion counts in a gene. Otherwise, the scatterplot is shown for individual TA sites
                -log           := If set, this shows the axes of the scatterplot on log scale
    """.replace("\n        ","\n     ")
    
    # usage removed above: {console_tools.subcommand_prefix} {cli_name} <combined_wig_file> <metadata_file> <annotation_file> --cond <comma-separated condition names> <output.png> [-genes -log]
    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=4)
        
        combined_wig_path      = args[0]
        metadata_path          = args[1]
        annotation_path        = args[2]
        output_path            = args[3] # png file
        avg_by_conditions = "cond" in kwargs
        condition_names = console_tools.string_arg_to_list(kwargs["cond"])
        sample_ids      = console_tools.string_arg_to_list(kwargs["samp"])
        
        combined_wig = tnseq_tools.CombinedWig.load(
            main_path=combined_wig_path,
            metadata_path=metadata_path,
            annotation_path=annotation_path,
        )
        # 
        # filter either by condition or by sample
        # 
        if avg_by_conditions:
            combined_wig = combined_wig.with_only(condition_names=condition_names)
        else:
            combined_wig = combined_wig.with_only(wig_ids=sample_ids)
        
        # save the data
        Method.output(
            combined_wig=combined_wig,
            annotation_path=annotation_path,
            output_path=output_path,
            avg_by_conditions=avg_by_conditions,
            normalization=kwargs["n"],
            gene_means="genes" in kwargs, # bool
            log_scale="log" in kwargs # bool
        )
        
    @gui.add_menu("Pre-Processing", "Visualize", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers, parameter_panel
        Method.value_getters = LazyDict()
        Method.by_condition = False
        with panel_helpers.NewPanel() as (panel, main_sizer):
            parameter_panel.set_instructions(
                title_text= self.name,
                sub_text= "",
                method_specific_instructions="""
                    A useful tool to show a detailed correlation of counts between 2 datasets.

                    1. Ensure the correct annotation file has been loaded in 

                    2. Select samples from the sample pane
                    
                    3. Select whether you would like to calculate the mean insertion count within a gene prior to the plot generation

                    4. [Optional] Select if you would like to have the plots be in log scale and/or normalized using a specific method

                    5. Click Run
                """.replace("\n                    ","\n")
            )

            sample_ids = [x.id for x in gui.samples]
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=Method.from_gui)
            Method.value_getters.combined_wig  = panel_helpers.combined_wig_filtered_by_sample_input(panel, main_sizer)
            Method.value_getters.gene_means    = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="average counts at the gene level", default_value=False, tooltip_text="if false, this shows the scatterplot of insertion counts at individual TA sites", widget_size=None)
            Method.value_getters.log_scale     = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="show axes on log scale", default_value=False, tooltip_text="show axes on log scale", widget_size=None)
            Method.value_getters.normalization = panel_helpers.create_normalization_input(panel, main_sizer)

    
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
            assert len(arguments.combined_wig.samples) >= 2, "Please select two or more samples on the left"
        else:
            assert len(arguments.combined_wig.condition_names) >= 2, "Please select only two or more conditions on the left"
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
            combined_wig = tnseq_tools.CombinedWig.load(main_path=combined_wig_path, metadata_path=metadata_path, annotation_path=annotation_path)
        
        with transit_tools.TimerAndOutputs(method_name=Method.identifier, output_paths=[output_path], disable=disable_logging) as timer:
            # 
            # by gene or site
            # 
            if gene_means:
                means, genes, labels = calc_gene_means(combined_wig=combined_wig, normalization=normalization, avg_by_conditions=avg_by_conditions, n_terminus=n_terminus, c_terminus=c_terminus)
                counts = means
            else:
                if combined_wig.metadata:
                    labels = combined_wig.metadata.wig_ids
                else:
                    labels = combined_wig.wig_fingerprints
                # average conditions if needed
                if avg_by_conditions:
                    combined_wig = combined_wig.averaged(by_conditions=True)
                    labels = combined_wig.metadata.condition_names
                counts = combined_wig.read_counts_array
            
            # 
            # plot
            # 
            import seaborn as sns
            import matplotlib.pyplot as plt
            import pandas as pd
            
            counts_df = pd.DataFrame(counts, columns=labels)
            plt.figure()
            if log_scale:  g = sns.PairGrid(numpy.log10(counts_df))
            else: g = sns.PairGrid(counts_df)
            g.map(sns.scatterplot, s=10, alpha=.5)
            
            if gene_means: g.fig.suptitle("scatter plot of mean insertion counts for each gene")
            else: g.fig.suptitle("scatter plot of insertion counts at individual TA sites")
            plt.savefig(output_path, bbox_inches='tight', pad_inches = 0.5)