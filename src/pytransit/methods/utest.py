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
from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.specific_tools.transit_tools import wx, basename


from pytransit.methods.pathway_enrichment import Method as PathwayEnrichment


@misc.singleton
class Method:
    name = "UTest"
    identifier  = name
    cli_name    = name.lower()
    description = f"""Mann-Whitney U-test of conditional essentiality"""
    menu_name   = f"{name} - {description}"
    himar1      = True # default assume
    
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
        "--n",
        "--iN",
        "--iC",
    ]
    usage_string = f"""
        Usage:
            {console_tools.subcommand_prefix} {cli_name} <combined_wig_file> <metadata_file> <annotation_file> <condition_for_control> <condition_for_experimental> <output_file> [Optional Arguments]

        Optional Arguments:
            --n <string>     :=  Normalization method. Default: --n TTR
            --iN <float>     :=  Ignore TAs occurring at given fraction (as integer) of the N terminus. Default: --iN 0
            --iC <float>     :=  Ignore TAs occurring at given fraction (as integer) of the C terminus. Default: --iC 0
    """.replace("\n        ", "\n")
    
    @gui.add_menu("Method", "Himar1", menu_name)
    def on_menu_click(event):
        Method.himar1 = True
        Method.define_panel(event)
    
    #@gui.add_menu("Method", "tn5", menu_name)
    #def on_menu_click(event):
    #    Method.himar1 = False
    #    Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers, parameter_panel
        with panel_helpers.NewPanel() as (panel, main_sizer):
            parameter_panel.set_instructions(
                title_text=self.name,
                sub_text="",
                method_specific_instructions="""
                    This is a method for comparing datasets from a TnSeq library evaluated in two different conditions, analogous to resampling. This is a rank-based test on whether the level of insertions in a gene or chromosomal region are significantly higher or lower in one condition than the other. Effectively, the insertion counts at the TA sites in the region are pooled and sorted. Then the combined ranks of the counts in region A are compared to those in region B, and p-value is calculated that reflects whether there is a significant difference in the ranks. The advantage of this method is that it is less sensitive to outliers (a unusually high insertion count at just a single TA site). A reference for this method is (Santa Maria et al., 2014).
                """.replace("\n                    ","\n"),
            )
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            self.value_getters = LazyDict()
            self.value_getters.control_condition        = panel_helpers.create_control_condition_input(panel, main_sizer)
            self.value_getters.experimental_condition   = panel_helpers.create_experimental_condition_input(panel, main_sizer)
            self.value_getters.n_terminus               = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus               = panel_helpers.create_c_terminus_input(panel, main_sizer)
            self.value_getters.normalization            = panel_helpers.create_normalization_input(panel, main_sizer)
            
    
    @staticmethod        
    def from_gui(frame):
        arguments = LazyDict()
        
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
            default_file_name=f"{Method.cli_name}_output.tsv",
            output_extensions='Common output extensions (*.tsv,*.dat,*.out)|*.txt;*.tsv;*.dat;*.out;|\nAll files (*.*)|*.*',
        )
        # if user didn't select an output path
        if not arguments.output_path:
            return None
        
        if arguments.control_condition == '[None]':
            logging.error(f'''Please select control condition''')
        if arguments.experimental_condition == '[None]':
            logging.error(f'''Please select experimental condition''')
        
        arguments.combined_wig = gui.combined_wigs[-1].with_only(condition_names=[ arguments.control_condition, arguments.experimental_condition ])
        Method.output(**arguments)

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=6)
        
        Method.output(
            combined_wig=tnseq_tools.CombinedWig.load(
                main_path=args[0],
                metadata_path=args[1],
                annotation_path=args[2],
            ),
            control_condition=args[3],
            experimental_condition=args[4],
            output_path=args[5],
            normalization=kwargs["n"],
            n_terminus=kwargs["iN"],
            c_terminus=kwargs["iC"],
        )
    
    @staticmethod
    def output(*, combined_wig, control_condition, experimental_condition, output_path, normalization=None, n_terminus=None, c_terminus=None, ignore_codon=None, significance_threshold=None, disable_logging=False):
        import scipy.stats
        from pytransit.specific_tools import stat_tools
        # Defaults (even if argument directly provided as None)
        normalization          = normalization          if normalization          is not None else "TTR"
        n_terminus             = n_terminus             if n_terminus             is not None else 0.0
        c_terminus             = c_terminus             if c_terminus             is not None else 0.0
        ignore_codon           = ignore_codon           if ignore_codon           is not None else True
        significance_threshold = significance_threshold if significance_threshold is not None else 0.05
        
        with transit_tools.TimerAndOutputs(method_name=Method.identifier, output_paths=[output_path], disable=disable_logging) as timer:
            number_of_control_wigs = sum( 1 for each in combined_wig.samples if control_condition in each.condition_names ) # TODO: check if this is right
            # 
            # restrict to relevent data
            # 
            combined_wig = combined_wig.with_only(condition_names=[ control_condition, experimental_condition ])
            
            # 
            # normalize
            # 
            if normalization != "nonorm":
                logging.log(f"Normalizing with {normalization}")
                combined_wig = combined_wig.normalized_with(normalization)
            
            control_samples_by_gene = combined_wig.with_only(condition_names=[control_condition]).get_genes(
                ignore_codon=ignore_codon,
                n_terminus=n_terminus,
                c_terminus=c_terminus,
            )
            experimental_samples_by_gene = combined_wig.with_only(condition_names=[experimental_condition]).get_genes(
                ignore_codon=ignore_codon,
                n_terminus=n_terminus,
                c_terminus=c_terminus,
            )

            # u-test
            data = []
            N = len(control_samples_by_gene)
            count = 0
            
            for progress, (control_gene, experimental_gene) in informative_iterator.ProgressBar(
                tuple(zip(control_samples_by_gene, experimental_samples_by_gene)),
                title=f"Running {Method.name}",
                disable_logging=True
            ):
                count += 1
                if (control_gene.k + experimental_gene.k) == 0 or control_gene.n == 0:
                    test_obs   = 0
                    mean1      = 0
                    mean2      = 0
                    log2_fc    = 0
                    u_stat     = 0.0
                    pval_2tail = 1.0
                else:
                    data1 = control_gene.reads.flatten()
                    data2 = experimental_gene.reads.flatten()
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
                        control_gene.orf,
                        control_gene.name,
                        control_gene.desc,
                        control_gene.n,
                        mean1,
                        mean2,
                        log2_fc,
                        u_stat,
                        pval_2tail,
                    ]
                )

                # # Update Progress
                # percent = (100.0 * count / N)
                # if gui.is_active:
                #     text = "Running Mann-Whitney U-test Method... %1.1f%%" % percent
                #     parameter_panel.progress_update(text, percent)

            data.sort()
            logging.log("")  # Printing empty line to flush stdout
            logging.log("Performing Benjamini-Hochberg Correction")
            qvals = stat_tools.bh_fdr_correction([row[Method.column_names.index("P Value")] for row in data])
            
            number_of_significant_genes = len([ 1 for each in qvals if each < significance_threshold ])
            
            # 
            # write output
            # 
            transit_tools.write_result(
                path=output_path, # path=None means write to STDOUT
                file_kind=Method.identifier,
                rows=[
                    [*row, qval]
                        for row, qval in zip(data, qvals) 
                ],
                column_names=Method.column_names,
                extra_info=dict(
                    calculation_time=f"{(timer.duration_in_seconds):0.1f}seconds",
                    analysis_type=Method.identifier,
                    stats=dict(
                        number_of_significant_genes=number_of_significant_genes,
                    ),
                    parameters={
                        "normalization": normalization,
                        "control_condition": control_condition,
                        "experimental_condition": experimental_condition,
                        "n_terminus":n_terminus,
                        "c_terminus":c_terminus,
                        "ignore_codon":ignore_codon,
                        "significance_threshold":significance_threshold,
                    },
                ),
            )


@transit_tools.ResultsFile
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier)
    
    def __init__(self, path=None):
        from pytransit.components.spreadsheet import SpreadSheet
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
                "Pathway Enrichment": lambda *args: PathwayEnrichment.call_from_results_panel(path),
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
            "summary": f"{number_of_significant} significant conditionally essential genes",
        })
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

