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
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    name = "Corrplot"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"Correlation Plot"
    description = f"""Make correlation plot"""
    
    transposons = [ "himar1" ] # not sure if this is right -- Jeff

    prev_menu_choice = None
    is_post_processing = False
    defaults = LazyDict(
        top_k=-1,
        q_value_threshold=0.05,
        low_mean_filter=5,
    )

    inputs = LazyDict(
        combined_wig=None,
        annotation_path=None,
        output_path=None,
        avg_by_conditions=False,
        normalization="TTR", #TRI hard-coded for now
    )
    
    valid_cli_flags = [
        "-avg_by_conditions",
        "--n",
        "--iN",
        "--iC",
        #if on output
        "-anova",
        "-zinb",
        "--qval",
        "--topk",
        "--low-mean-filter",
        
    ]
    #TRI - could add a flag for Adj P Value cutoff (or top n most signif genes)

    usage_string = f"""
        Usage:
        As a Pre-processing visualization:

            {console_tools.subcommand_prefix} {cli_name} <combined_wig_file> <metadata_file> <annotation_file> <output.png> [Optional Arguments]
        Optional Arguments:
            -avg_by_conditions := groups by conditions, take the mean, then show correlation between conditions. Default: false
            --n <string> := Normalization method. Default: --n TTR
            --iN    <N>  := Ignore TAs within given percentage (e.g. 5) of N terminus. Default: --iN 0
            --iC    <N>  := Ignore TAs within given percentage (e.g. 5) of C terminus. Default: --iC 0

        
        As a post processing visualization on ANOVA and ZINB output (avg_by_conditions = True)
        Usage 1:
            {console_tools.subcommand_prefix} {cli_name} <combined_wig_file> <metadata_file> <annotation_file> <anova_output> <output.png> -anova [Optional Arguments]
        Usage 2:
            {console_tools.subcommand_prefix} {cli_name}<combined_wig_file> <metadata_file> <annotation_file>  <zinb_output> <heatmap.png> -zinb [Optional Arguments]
        
        Optional Arguments:
            --topk <int>            := number of results
            --qval <float>          := adjusted p value threshold. Default --qval 0.05
            --low-mean-filter <int> := Filter out genes with grand mean count (across all conditions) below this threshold
                                    (even if adjusted p-value < 0.05)
                                    Default --low-mean-filter 5
    """
    
    # 
    # CLI method
    # 
    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=3)


        if not kwargs["anova"] and not kwargs["zinb"]:
            # map data to the core function
            Method.output(
                combined_wig=tnseq_tools.CombinedWig.load(
                    main_path=args[0],
                    metadata_path=args[1],
                    annotation_path=args[2],
                ),
                normalization=kwargs["n"],
                n_terminus=kwargs["iN"],
                c_terminus=kwargs["iC"],
                avg_by_conditions="avg_by_conditions" in kwargs, # bool
                output_path=args[3],
                disable_logging=False,
            )
        else:
            Method.load_from(
                filetype="anova" if kwargs["anova"] else "zinb",
                combined_wig_path=args[0],
                metadata_path=args[1],
                annotation_path=args[2],
                input_path        = args[3],
                output_path       = args[4],
                q_value_threshold = kwargs["qval"],
                top_k             = kwargs["topk"],
                low_mean_filter   = kwargs["low-mean-filter"], # filter out genes with grandmean<5 by default
            )

    @staticmethod
    def load_from(filetype, combined_wig_path,metadata_path, annotation_path, input_path, output_path, top_k=None, q_value_threshold=None, low_mean_filter=None):
        if filetype == "anova":
            from pytransit.methods.anova import File as AnovaFile
            AnovaFile(path=input_path).create_corrplot(output_path, combined_wig_path,metadata_path, annotation_path, top_k=top_k, qval=q_value_threshold, low_mean_filter=low_mean_filter)
        elif filetype == "zinb":
            from pytransit.methods.zinb import File as ZinbFile
            ZinbFile(path=input_path).create_corrplot(output_path, combined_wig_path,metadata_path, annotation_path, topk=top_k, qval=q_value_threshold, low_mean_filter=low_mean_filter)
      
    
    # 
    # Panel method
    # 
    @gui.add_menu("Pre-Processing", "Visualize", menu_name)
    def on_menu_click(event):
        Method.define_pre_processing_panel(event)
    
    def define_pre_processing_panel(self, _):
        from pytransit.components import panel_helpers, parameter_panel
        with panel_helpers.NewPanel() as (panel, main_sizer):
            parameter_panel.set_instructions(
                title_text= self.name,
                sub_text= "",
                method_specific_instructions="""
                    A useful tool when evaluating the quality of a collection of TnSeq datasets is to make a correlation plot of the mean insertion counts (averaged at the gene-level) among samples.

                    1. Ensure the correct annotation file has been loaded in 
                    
                    2. Select whether you would like to calculate the means across replicates within a condition

                    3. Click Run
                """.replace("\n                    ","\n")
            )

            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            self.value_getters = LazyDict()
            self.value_getters.avg_by_conditions = panel_helpers.create_check_box_getter(panel, main_sizer, label_text="average counts by condition", default_value=False, tooltip_text="correlations among conditions (where counts are averaged among replicates of each condition) versus all individual samples", widget_size=None)
            self.value_getters.normalization     = panel_helpers.create_normalization_input(panel, main_sizer)

    @gui.add_menu("Post-Processing", "ANOVA", "Corrplot")
    def on_menu_click(event):
        Method.prev_menu_choice = "anova"
        Method.define_postprocessing_panel(event)
    
    @gui.add_menu("Post-Processing", "ZINB", "Corrplot")
    def on_menu_click(event):
        Method.prev_menu_choice = "zinb"
        Method.define_postprocessing_panel(event)

    def define_postprocessing_panel(self, _):
        Method.is_post_processing = True
        from pytransit.components import panel_helpers, parameter_panel
        with panel_helpers.NewPanel() as (panel, main_sizer):
            parameter_panel.set_instructions(
                title_text=self.name,
                sub_text="",
                method_specific_instructions="""
                    The output of ANOVA or ZINB can be used to generate a corrplot that shows the correlation of the original counts of significant genes conditions-wise, which is s useful for shedding light on the relationships among the conditions apparent in the data.
                """.replace("\n                    ","\n"),
            )
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            
            self.value_getters = LazyDict()
            self.value_getters.input_path        = panel_helpers.create_file_input(  panel, main_sizer, button_label=f"Select {Method.prev_menu_choice} file", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*')
            self.value_getters.q_value_threshold = panel_helpers.create_float_getter(panel, main_sizer, label_text="Adj P Value Cutoff", default_value=Method.defaults.q_value_threshold, tooltip_text="Change adjusted p-value threshold for selecting genes")
            self.value_getters.top_k             = panel_helpers.create_int_getter(  panel, main_sizer, label_text="Top K",              default_value=Method.defaults.top_k,             tooltip_text="(-1 means all) Sometimes there are so many genes it is hard to see the heatmap top genes. This allows limiting to the top K genes (ranked by significance; adjusted p-value)")
            self.value_getters.low_mean_filter   = panel_helpers.create_float_getter(panel, main_sizer, label_text="Low Mean Filter",    default_value=Method.defaults.low_mean_filter,   tooltip_text="Filter out genes with grand mean count (across all conditions) below this threshold (even if adjusted p-value < 0.05)")
            self.value_getters.normalization     = panel_helpers.create_normalization_input(panel, main_sizer)
            self.value_getters.n_terminus      = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus      = panel_helpers.create_c_terminus_input(panel, main_sizer)

    @staticmethod
    def from_gui(frame):
        combined_wig = gui.combined_wigs[-1]
        arguments = LazyDict()
        # 
        # get global data
        # 
        arguments.combined_wig = gui.combined_wigs[-1] #TRI what if not defined? fail gracefully?
        
        # 
        # call all GUI getters, puts results into respective Method.defaults key-value
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
            default_file_name=f"corrplot.png",
            output_extensions='PNG file (*.png)|*.png;|\nAll files (*.*)|*.*',
        )
        
        # if user didn't select an output path
        if not arguments.output_path:
            return None
        
        if Method.is_post_processing:
            Method.load_from(
                filetype          = Method.prev_menu_choice,
                combined_wig_path = combined_wig.main_path,
                metadata_path     = combined_wig.metadata_path,
                annotation_path   = gui.annotation_path,
                input_path        = arguments.input_path,
                output_path       = arguments.output_path,
                q_value_threshold = arguments.q_value_threshold,
                top_k             = arguments.top_k,
                low_mean_filter   = arguments.low_mean_filter, # filter out genes with grandmean<5 by default
            )
        else:
        # run the core function directly
            Method.output(
                combined_wig_path = combined_wig.main_path,
                metadata_path     = combined_wig.metadata_path,
                annotation_path   = gui.annotation_path,
                avg_by_conditions=arguments.avg_by_conditions, # bool
                output_path=arguments.output_path,
                disable_logging=False,
            )
    
    corrplot_r_function = None    
    @staticmethod
    def output(*, combined_wig_path=None, metadata_path=None, annotation_path=None, combined_wig=None, normalization=None, formatted_rows=None, top_k=None, q_value_threshold=None, low_mean_filter=None, avg_by_conditions=None, output_path=None, n_terminus=None, c_terminus=None, disable_logging=False):
        # Defaults (even if argument directly provided as None)
        normalization     = normalization     if normalization     is not None else "TTR"
        avg_by_conditions = avg_by_conditions if avg_by_conditions is not None else False
        output_path       = output_path       if output_path       is not None else "./corrplot.png"
        n_terminus        = n_terminus        if n_terminus        is not None else 0.0
        c_terminus        = c_terminus        if c_terminus        is not None else 0.0
        top_k              = int(top_k)               if top_k              != None  else Method.defaults.top_k
        q_value_threshold  = float(q_value_threshold) if q_value_threshold  != None  else Method.defaults.q_value_threshold
        low_mean_filter    = int(low_mean_filter)     if low_mean_filter    != None  else Method.defaults.low_mean_filter
        
        if combined_wig == None:
            combined_wig = tnseq_tools.CombinedWig.load(main_path=combined_wig_path,metadata_path=metadata_path, annotation_path=annotation_path)
        
        from pytransit.methods.gene_means import Method as GeneMeansMethod
        with transit_tools.TimerAndOutputs(method_name=Method.identifier, output_paths=[output_path], disable=disable_logging,):
            import seaborn as sns
            import matplotlib.pyplot as plt
            import numpy as np
            import pandas as pd

            _, (means, genes, headers) = GeneMeansMethod.calculate(combined_wig, normalization=normalization, avg_by_conditions=avg_by_conditions, n_terminus=n_terminus, c_terminus=c_terminus)
            
            position_hash = {}
            for i, col in enumerate(headers):
                position_hash[col] = [x[i] for x in means]
            df = pd.DataFrame.from_dict(position_hash, orient="columns") 
            df.index=[genes[i].orf+"/"+genes[i].name for i in range(len(df))]

            clean_columns = [s.split("/")[-1].split(".wig")[0] for s in df.columns]
            clean_headers = [s.split("/")[-1].split(".wig")[0] for s in headers]
            df.columns = clean_columns

            logging.log("Filtering Rows to include")
            if formatted_rows != None:
                sorted_rows = sorted(formatted_rows, key=lambda row: row["q_value"]) 
                slice_end = len(sorted_rows) if top_k == -1 else top_k
                sorted_rows = sorted_rows[:slice_end]
                significant_rows = [ each for each in sorted_rows if each["q_value"] < q_value_threshold ]
                if len(significant_rows) >= top_k: sorted_rows = significant_rows
                gene_names = []
                for each_row in sorted_rows:
                    mean_of_means = round(numpy.mean(each_row["means"]), 1)
                    if mean_of_means < low_mean_filter:
                        print(f"""excluding {each_row["gene_name"]}, mean(means)={mean_of_means}""")
                    else:
                        gene_names.append(each_row["gene_name"])

                df = df[df.index.isin(gene_names)]

            logging.log("Creating the Correlation Plot")
            corr_df = df[clean_headers].corr()
            lower_corr_df = corr_df.where(np.tril(np.ones(corr_df.shape)).astype(np.bool)) #only plot the bottom triangle
            

            plt.figure()
            g = sns.heatmap(lower_corr_df, cmap=sns.color_palette('blend:Red,White,Blue', as_cmap=True),  square=True, linewidths=.5, cbar_kws={"shrink": .5}, vmin=-1, vmax=1)
            if Method.is_post_processing:
                plt.title("Correlations of Genes in Conditions in output of "+Method.prev_menu_choice)
            elif avg_by_conditions == False:
                plt.title("Correlations of Genes in Samples")
            else:
                plt.title("Correlations of Genes in Conditions")

            if os.path.exists(output_path): os.remove(output_path)
            g.figure.savefig(output_path, bbox_inches='tight')