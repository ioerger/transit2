import sys
import os
import time
import ntpath
import math
import random
import datetime
import heapq
import collections

import numpy
from pytransit.generic_tools.lazy_dict import LazyDict

from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.components.parameter_panel import panel as parameter_panel
from pytransit.components.parameter_panel import progress_update
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.specific_tools import  gui_tools, transit_tools, console_tools, tnseq_tools, norm_tools
from pytransit.generic_tools import csv, misc
import pytransit.components.results_area as results_area

@misc.singleton
class Method:
    identifier  = "Tnseq Stats"
    cli_name    = identifier.replace(" ","_").lower()
    identifier  = "TnseqStats"
    menu_name   = f"Summary statistics ({identifier})"
    description = """Analyze statistics of TnSeq datasets in combined_wig file"""
    
    inputs = LazyDict(
        combined_wig=None,
        normalization="nonorm",
        output_path=None,
    )
    
    column_names = [
        "Dataset",
        "Density",
        "Mean Count",
        "Non Zero Mean",
        "Non Zero Median",
        "Max Count",
        "Total Counts",
        "Skewness",
        "Kurtosis",
        "Pickands Tail Index",
    ]
    
    valid_cli_flags = [
        "--n", # normalization flag
    ]
    
    usage_string = f"""
        Usage: {console_tools.subcommand_prefix} tnseq_stats <wig_file or combined_wig_file> [<output_file>]
    """
    
    @gui.add_menu("Pre-Processing", "Generate", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            # only need Norm selection and Run button        
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            self.value_getters = LazyDict(
                normalization=panel_helpers.create_normalization_input(panel, main_sizer,default="nonorm")
            )

    @staticmethod
    def from_gui(frame):
        # 
        # get wig files
        # 
        combined_wig = gui.combined_wigs[-1]
        Method.inputs.combined_wig = combined_wig.main_path
        
        # 
        # setup custom inputs
        # 
        for each_key, each_getter in Method.value_getters.items():
            try:
                Method.inputs[each_key] = each_getter()
            except Exception as error:
                raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')

        # 
        # save result files
        # 
        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name="tnseq_stats.tsv",
            output_extensions='Common output extensions (*.tsv,*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
        )
        if not Method.inputs.output_path:
            return None

        return Method

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        if len(args)==0: print(Method.usage_string); sys.exit() # TRI
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        
        wig_or_combined_wig = console_tools.string_arg_to_list(args[0])
        output_path = args[1] if len(args) > 1 else None
        normalization = kwargs.get("n", Method.inputs.normalization)

        if len(wig_or_combined_wig) > 1:
            wigs = wig_or_combined_wig
            combined_wig = None
        else:
            from pytransit.methods.combined_wig import Method as CombinedWigMethod
            
            wig_or_combined_wig = wig_or_combined_wig[0]
            if CombinedWigMethod.file_is_combined_wig(wig_or_combined_wig):
                combined_wig = wig_or_combined_wig
                wigs = []
            else:
                combined_wig = None
                wigs = [ wig_or_combined_wig ]
        
        Method.inputs.update(dict(
            wigs=wigs,
            combined_wig=combined_wig, 
            normalization=normalization,
            output_path=output_path,
        ))
        
        Method.Run()
        
    def Run(self):
        logging.log("Starting tnseq_stats analysis")
        start_time = time.time()

        # 
        # get data
        # 
        logging.log(f"Getting Data from {self.inputs.combined_wig}")
        if self.inputs.combined_wig:
            combined_wig = tnseq_tools.CombinedWig.load(main_path=self.inputs.combined_wig)
            ta_site_positions, read_counts, filenames_in_comb_wig = combined_wig.as_tuple
        else:
            read_counts, ta_site_positions = tnseq_tools.CombinedWig.gather_wig_data(self.inputs.wigs)
            filenames_in_comb_wig = self.inputs.wigs
            
        logging.log(f"Normalizing using: {self.inputs.normalization}")
        read_counts, factors = norm_tools.normalize_data(read_counts, self.inputs.normalization)
            
        # 
        # process read_counts
        # 
        logging.log("processing read_counts")
        results = self.calc_tnseq_stats(read_counts, filenames_in_comb_wig)

        # 
        # write output
        # 
        if True:
            # 
            # write to file
            # 
            transit_tools.write_result(
                path=self.inputs.output_path, # path=None means write to STDOUT
                file_kind=Method.identifier,
                rows=[ list(row) for row in results ],
                column_names=self.column_names,
                extra_info=dict(
                    parameters=dict(
                        normalization=self.inputs.normalization,
                    ),
                ),
            )
            if self.inputs.output_path != None:
                logging.log(f"Adding File: {self.inputs.output_path}")
                results_area.add(self.inputs.output_path)
        
        logging.log("Finished TnseqStats")
        logging.log("Time: %0.1fs\n" % (time.time() - start_time))

    def pickands_tail_index(self, vals):
        srt = sorted(vals, reverse=True)
        PTIs = []
        for M in range(10, 100):
            PTI = numpy.log(
                (srt[M] - srt[2 * M]) / float(srt[2 * M] - srt[4 * M])
            ) / numpy.log(2.0)
            PTIs.append(PTI)
        return numpy.median(PTIs)

    # data=numpy array of (normalized) insertion counts at TA ta_site_positions for multiple samples; sample_names=.wig filenames
    def calc_tnseq_stats(self, read_counts, sample_names): 
        results = []
        for i in range(read_counts.shape[0]):
            (
                density,
                meanrd,
                nzmeanrd,
                nzmedianrd,
                maxrd,
                totalrd,
                skew,
                kurtosis,
            ) = tnseq_tools.get_data_stats(read_counts[i, :])
            nzmedianrd = int(nzmedianrd) if numpy.isnan(nzmedianrd) == False else 0
            pti = self.pickands_tail_index(read_counts[i, :])
            vals = [
                sample_names[i],
                "%0.3f" % density,
                "%0.1f" % meanrd,
                "%0.1f" % nzmeanrd,
                "%d" % nzmedianrd,
                maxrd,
                int(totalrd),
                "%0.1f" % skew,
                "%0.1f" % kurtosis,
                "%0.3f" % pti
            ]
            results.append(vals)
        return results


@transit_tools.ResultsFile
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=transit_tools.basename(self.path),
            type=Method.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title="tnseq_stats",heading="",column_names=self.column_names,rows=self.rows).Show(),
            })
        )
        
        # 
        # get column names
        # 
        comments, headers, rows = csv.read(self.path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
        if len(comments) == 0:
            raise Exception(f'''No comments in file, and I expected the last comment to be the column names, while to load tnseq_stats file "{self.path}"''')
        self.column_names = comments[-1].split("\t")
        
        # 
        # get rows
        #
        self.rows = []
        for each_row in rows:
            row = {}
            for each_column_name, each_cell in zip(self.column_names, each_row):
               row[each_column_name] = each_cell
            self.rows.append(row)
        
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
    
    
