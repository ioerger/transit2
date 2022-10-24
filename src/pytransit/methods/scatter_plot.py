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
    cli_name = name.replace(" ",'').lower()
    
    # save the data
    inputs = dict(
        avg_by_conditions=False,
        combined_wig=None,
        annotation_path=None,
        output_path=None,
        normalization="TTR",
        gene_means=False,
        log_scale=False,
        n_terminus=0,
        c_terminus=0,
    )
    
    valid_cli_flags = [ "log", "genes", "cond" ]

    usage_string = f"""
        usage:
            {console_tools.subcommand_prefix} {cli_name} <combined_wig> <metadata_file> <sample_id1> <sample_id2> <annotation_file> <output.png> [-genes -log]
            {console_tools.subcommand_prefix} {cli_name} <combined_wig> <metadata_file> -cond <condition1> <condition2> <annotation_file> <output.png> [-genes -log]
    """.replace("\n        ","\n     ")
    
    @gui.add_menu("Pre-Processing", menu_name+" By Sample")
    def on_menu_click(event):
        from pytransit.components import panel_helpers
        self.value_getters = LazyDict()
        self.inputs.avg_by_conditions = False
        with panel_helpers.NewPanel() as (panel, main_sizer):
            sample_ids = [x.id for x in gui.samples]
            self.value_getters.combined_wig  = panel_helpers.combined_wig_filtered_by_sample_input(panel, main_sizer),
            self.value_getters.gene_means    = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="average counts at the gene level", default_value=False, tooltip_text="if false, this shows the scatterplot of insertion counts at individual TA sites", widget_size=None)
            self.value_getters.log_scale     = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="show axes on log scale", default_value=False, tooltip_text="show axes on log scale", widget_size=None)
            self.value_getters.normalization = panel_helpers.create_normalization_input(panel, main_sizer)
            self.value_getters.n_terminus    = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus    = panel_helpers.create_c_terminus_input(panel, main_sizer)

            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
    
    @gui.add_menu("Pre-Processing", menu_name+" By Condition")
    def on_menu_click(event):
        from pytransit.components import panel_helpers
        self.value_getters = LazyDict()
        self.inputs.avg_by_conditions = True
        with panel_helpers.NewPanel() as (panel, main_sizer):
            sample_ids = [x.id for x in gui.samples]
            self.value_getters.combined_wig  = panel_helpers.combined_wig_filtered_by_condition_input(panel, main_sizer),
            self.value_getters.gene_means    = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="average counts at the gene level", default_value=False, tooltip_text="if false, this shows the scatterplot of insertion counts at individual TA sites", widget_size=None)
            self.value_getters.log_scale     = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="show axes on log scale", default_value=False, tooltip_text="show axes on log scale", widget_size=None)
            self.value_getters.normalization = panel_helpers.create_normalization_input(panel, main_sizer)
            self.value_getters.n_terminus    = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus    = panel_helpers.create_c_terminus_input(panel, main_sizer)

            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
    
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

        #        # determine which 2 samples were selected in samples area...
        #        #TRI if 3 or more samples selected, show grid of scatterplots?
        #        if (len(gui.selected_samples))<2: 
        #          logging.error("need to select at least 2 samples for making scatter plot")
        #        metadata = tnseq_tools.CombinedWigMetadata(Method.inputs.metadata)
        #        fp1 = gui.selected_samples[0].fingerprint
        #        Method.inputs.sample1 = metadata.id_for(fp1)
        #        fp2 = gui.selected_samples[1].fingerprint
        #        Method.inputs.sample2 = metadata.id_for(fp2)
        
        # if plotting samples, (for now) only allow two samples
        if not self.inputs.avg_by_conditions:
            assert len(self.inputs.combined_wig.samples) == 2, "Please select only two samples on the left"
        else:
            assert len(self.inputs.combined_wig.condition_names) == 2, "Please select only two conditions on the left"
        
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
        
        combined_wig_path      = args[0]
        metadata_path          = args[1]
        sample_or_condition1   = args[2]
        sample_or_condition2   = args[3]
        annotation_path        = args[4]
        output_path            = args[5] # png file
        avg_by_conditions = "cond" in kwargs
        
        combined_wig = tnseq_tools.CombinedWig(
            main_path=combined_wig,
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
        Method.inputs.update(dict(
            combined_wig = combined_wig,
            annotation_path = annotation_path,
            output_path = output_path,
            normalization = kwargs.get("n", Method.inputs.normalization),
            gene_means = "genes" in kwargs, # bool
            log_scale = "log" in kwargs # bool
        ))
        
        Method.Run()
    
    def Run(self):
        logging.log(f"Starting {Method.identifier} analysis")
        start_time = time.time()
        
        transit_tools.make_scatterplot(
            combined_wig=self.inputs.combined_wig,
            normalization=self.inputs.normalization,
            annotation_path=self.inputs.annotation_path,
            avg_by_conditions=self.inputs.avg_by_conditions,
            output_path=self.inputs.output_path,
            gene_means = self.inputs.gene_means,
            n_terminus=self.inputs.n_terminus,
            c_terminus = self.inputs.c_terminus,
            log_scale = self.inputs.log_scale
        )
        
        if gui.is_active:
            logging.log(f"Adding File: {self.inputs.output_path}")
            results_area.add(self.inputs.output_path)
        logging.log(f"Finished {Method.identifier} analysis in {time.time() - start_time:0.1f}sec")