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
from pytransit.basics.lazy_dict import LazyDict

from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components.parameter_panel import panel as parameter_panel
from pytransit.components.parameter_panel import progress_update
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.tools import logging, gui_tools, transit_tools, console_tools, tnseq_tools, norm_tools
from pytransit.basics import csv, misc
import pytransit.components.results_area as results_area



@misc.singleton
class Analysis:
    identifier  = "tnseq_stats"
    cli_name    = identifier.lower()
    menu_name   = f"{identifier} - Analyze statistics of TnSeq datasets"
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
        "-n", # normalization flag
        "-o", # output filename (optional)
        "-c", # indicates whether input is list of wig files (comma- or space-separated?), or a combined_wig file
    ]
    usage_string = f"""usage: {console_tools.subcommand_prefix} tnseq_stats <file.wig>+ [-o <output_file>]\n       {console_tools.subcommand_prefix} tnseq_stats -c <combined_wig> [-o <output_file>]"""
    
    @gui.add_menu("Analysis", "himar1", menu_name)
    def on_menu_click(event):
        Analysis.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            # only need Norm selection and Run button        
            self.value_getters = LazyDict(
                normalization=panel_helpers.create_normalization_input(panel, main_sizer,default="nonorm")
            )
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)

    @staticmethod
    def from_gui(frame):
        # 
        # get wig files
        # 
        combined_wig = gui.combined_wigs[0]
        Analysis.inputs.combined_wig = combined_wig.main_path
        
        # 
        # setup custom inputs
        # 
        for each_key, each_getter in Analysis.value_getters.items():
            try:
                Analysis.inputs[each_key] = each_getter()
            except Exception as error:
                raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')

        # 
        # save result files
        # 
        Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name="tnseq_stats.dat",
            output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
        )
        if not Analysis.inputs.output_path:
            return None

        return Analysis

    @staticmethod
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Analysis.usage_string)
        console_tools.handle_unrecognized_flags(Analysis.valid_cli_flags, kwargs, Analysis.usage_string)

        wigs = args # should be args[0]?
        combined_wig  = kwargs.get("c", Analysis.inputs.combined_wig)
        normalization = kwargs.get("n", Analysis.inputs.normalization) 
        output_path   = kwargs.get("o", Analysis.inputs.output_path)

        if combined_wig == None and len(wigs) == 0:
            logging.error(Analysis.usage_string)

        # save all the data
        Analysis.inputs.update(dict(
            wigs=wigs, ### what if user gives a list of wig files instead of a combined_wig?
            combined_wig=combined_wig, 
            normalization=normalization,
            output_path=output_path,
        ))
        
        return Analysis
        
    def Run(self):
        logging.log("Starting tnseq_stats analysis")
        start_time = time.time()

        # if you want to see which samples were selected...
        #from pytransit.components.samples_area import sample_table
        #datasets_selected = [ each_row["path"] for each_row in sample_table.selected_rows ]
        #for x in datasets_selected: print(str(x))

        # 
        # get data
        # 
        logging.log(f"Getting Data from {self.inputs.combined_wig}")
        sites, data, filenames_in_comb_wig = tnseq_tools.read_combined_wig(self.inputs.combined_wig)
        logging.log(f"Normalizing using: {self.inputs.normalization}")
        data, factors = norm_tools.normalize_data(data, self.inputs.normalization)
            
        # 
        # process data
        # 
        logging.log("processing data")
        results = self.calc_tnseq_stats(data,filenames_in_comb_wig)

        # 
        # write output
        # 
        if True:
            # note: first comment line is filetype, last comment line is column headers
            file = sys.stdout # print to console if not output file defined
            if self.inputs.output_path != None:
                file = open(self.inputs.output_path, "w")
            file.write("#%s\n" % self.identifier)
            file.write("#normalization: %s\n" % self.inputs.normalization)
            file.write("#"+"\t".join(self.column_names)+"\n")

            for vals in results:
                file.write("\t".join([str(x) for x in vals]) + "\n")

            if self.inputs.output_path != None:
                file.close()
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

    # data=numpy array of (normalized) insertion counts at TA sites for multiple samples; sample_names=.wig filenames
    def calc_tnseq_stats(self,data,sample_names): 
        results = []
        for i in range(data.shape[0]):
            (
                density,
                meanrd,
                nzmeanrd,
                nzmedianrd,
                maxrd,
                totalrd,
                skew,
                kurtosis,
            ) = tnseq_tools.get_data_stats(data[i, :])
            nzmedianrd = int(nzmedianrd) if numpy.isnan(nzmedianrd) == False else 0
            pti = self.pickands_tail_index(data[i, :])
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
        return transit_tools.file_starts_with(path, '#'+Analysis.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=transit_tools.basename(self.path),
            type=Analysis.identifier,
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
            File for {Analysis.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
    
    
Method = GUI = Analysis # for compatibility with older code/methods