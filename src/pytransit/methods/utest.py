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
    name = "UTest"
    identifier  = name
    cli_name    = name.lower()
    description = f"""Mann-Whitney U-test of conditional essentiality"""
    menu_name   = f"{name} - {description}"
    
    
    inputs = LazyDict(
        ctrldata=None,
        expdata=None,
        annotation_path=None,
        output_path=None,
        normalization="TTR",
        include_zeros=False,
        LOESS=False,
        ignore_codon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        significance_threshold=0.05,
    )
    
    column_names = [
        "Orf",
        "Name",
        "Desc",
        "Sites",
        "Mean Ctrl",
        "Mean Exp",
        "Log 2 FC",
        "U Statistic",
        "P Value",
        "Adj P Value",
    ]
    
    valid_cli_flags = [
        "-n",
        "-iz",
        "-l",
        "-iN",
        "-iC",
    ]
    usage_string = f"""
        Usage: {console_tools.subcommand_prefix} {cli_name} <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]

        Optional Arguments:
        -n <string>     :=  Normalization method. Default: -n TTR
        -iz             :=  Include rows with zero accross conditions.
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring at given fraction (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring at given fraction (as integer) of the C terminus. Default: -iC 0
    """.replace("\n        ", "\n")
    
    @gui.add_menu("Method", "himar1", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    @gui.add_menu("Method", "tn5", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            parameter_panel.set_instructions(
                title_text=self.name,
                sub_text="",
                method_specific_instructions="""
                    HANDLE_THIS
                """.replace("\n                    ","\n"),
            )
            self.value_getters = LazyDict()
            # panel_helpers.create_float_getter(panel, main_sizer, label_text="", default_value=0, tooltip_text="")
            # panel_helpers.create_int_getter(panel, main_sizer, label_text="", default_value=0, tooltip_text="")
            # panel_helpers.create_file_input(panel, main_sizer, button_label="", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*')
            # panel_helpers.create_choice_input(panel, main_sizer, label="", options=[], default_option=None, tooltip_text="")
            # panel_helpers.create_text_box_getter(panel, main_sizer, label_text="", default_value="", tooltip_text="", label_size=None, widget_size=None,)
            # panel_helpers.create_check_box_getter(panel, main_sizer, label_text="", default_value=False, tooltip_text="", widget_size=None)
            # @panel_helpers.create_button(panel, main_sizer, label="")
            # def when_button_clicked(event):
            #     print("do stuff")
            # @panel_helpers.create_button(panel, main_sizer, label="Show pop up")
            # def when_button_clicked(event):
            #     from pytransit.components import pop_up
            #     @pop_up.create_pop_up(panel)
            #     def create_pop_up_contents(pop_up_panel, sizer, refresh, close):
            # 
            #         @panel_helpers.create_button(pop_up_panel, sizer, label="Click me for pop up")
            #         def when_button_clicked(event):
            #             print("do stuff")
            
            self.value_getters.n_terminus             = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus             = panel_helpers.create_c_terminus_input(panel, main_sizer)
            self.value_getters.normalization          = panel_helpers.create_normalization_input(panel, main_sizer)
            
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
        
        # 
        # ask for output path(s)
        # 
        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}_output.tsv",
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
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=4)
        
        # FIXME: create combined wig object with only relevent data
        
        # save the data
        Method.inputs.update(dict(
            ctrldata = args[0].split(","),
            expdata = args[1].split(","),
            annotation_path = args[2],
            output_path = args[3],
            normalization = kwargs.get("n", "TTR"),
            include_zeros = kwargs.get("iz", False),
            LOESS = kwargs.get("l", False),
            ignore_codon = True,
            n_terminus = float(kwargs.get("iN", 0.00)),
            c_terminus = float(kwargs.get("iC", 0.00)),
        ))
        
        Method.Run()
    
    def Run(self):
        import scipy.stats
        from pytransit.specific_tools import stat_tools
        logging.log("Starting Mann-Whitney U-test Method")
        start_time = time.time()
        
        number_of_control_wigs = 1 # TODO: check if this is right
        
        # 
        # normalize
        # 
        if self.inputs.normalization != "nonorm":
            logging.log(f"Normalizing with {self.inputs.normalization}")
            combined_wig = combined_wig.normalized_with(self.inputs.normalization)
        
        if self.inputs.LOESS:
            logging.log("Performing LOESS Correction")
            combined_wig = combined_wig.with_loess_correction()

        G = combined_wig.get_genes(
            ignore_codon=self.inputs.ignore_codon,
            n_terminus=self.inputs.n_terminus,
            c_terminus=self.inputs.c_terminus,
        )

        # u-test
        data = []
        N = len(G)
        count = 0
        
        for gene in G:
            count += 1
            if gene.k == 0 or gene.n == 0:
                (test_obs, mean1, mean2, log2_fc, u_stat, pval_2tail) = (
                    0,
                    0,
                    0,
                    0,
                    0.0,
                    1.00,
                )
            else:

                if not self.inputs.include_zeros:
                    ii = numpy.sum(gene.reads, 0) > 0
                else:
                    ii = numpy.ones(gene.n) == 1

                data1 = gene.reads[:number_of_control_wigs, ii].flatten()
                data2 = gene.reads[number_of_control_wigs:, ii].flatten()
                try:
                    u_stat, pval_2tail = scipy.stats.mannwhitneyu(
                        data1, data2, alternative="two-sided"
                    )
                except ValueError as e:
                    u_stat, pval_2tail = 0.0, 1.00

                n1 = len(data1)
                n2 = len(data2)

                mean1 = 0
                if n1 > 0:
                    mean1 = numpy.mean(data1)
                mean2 = 0
                if n2 > 0:
                    mean2 = numpy.mean(data2)

                try:
                    # Only adjust log2_fc if one of the means is zero
                    if mean1 > 0 and mean2 > 0:
                        log2_fc = math.log((mean2) / (mean1), 2)
                    else:
                        log2_fc = math.log((mean2 + 1.0) / (mean1 + 1.0), 2)
                except:
                    log2_fc = 0.0

            data.append(
                [
                    gene.orf,
                    gene.name,
                    gene.desc,
                    gene.n,
                    mean1,
                    mean2,
                    log2_fc,
                    u_stat,
                    pval_2tail,
                ]
            )

            # Update Progress
            percent = (100.0 * count / N)
            text = "Running Mann-Whitney U-test Method... %1.1f%%" % percent
            parameter_panel.progress_update(text, percent)

        logging.log("")  # Printing empty line to flush stdout
        logging.log("Performing Benjamini-Hochberg Correction")
        data.sort()
        qval = stat_tools.bh_fdr_correction([row[Method.column_names.index("P Value")] for row in data])
        
        number_of_significant_genes = len([ 1 for each in qval if each > self.inputs.significance_threshold ])
        
        # 
        # write output
        # 
        transit_tools.write_result(
            path=self.inputs.output_path, # path=None means write to STDOUT
            file_kind=Method.identifier,
            rows=[
                [*row, qval]
                    for row, qval in zip(data, qval) 
            ],
            column_names=Method.column_names,
            extra_info=dict(
                stats=dict(
                    number_of_significant_genes=number_of_significant_genes,
                ),
                parameters=self.inputs,
            ),
        )

        logging.log("Adding File: %s" % (self.inputs.output_path))
        results_area.add(self.inputs.output_path)
        logging.log("Finished Mann-Whitney U-test Method")

@transit_tools.ResultsFile
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier,
                    heading=self.comments_string or misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[
                        "Adj P Value",
                    ],
                ).Show(),
            })
        )
        
        # 
        # read in data
        # 
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        
        # 
        # get summary stats
        #
        number_of_significant = self.extra_data["stats"]["number_of_significant_genes"]
        self.values_for_result_table.update({
            " ": f"{number_of_significant} significant conditionally essential genes",
        })
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

