import sys
import os
import time
import numpy

from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import misc, informative_iterator
from pytransit.generic_tools import csv
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

@misc.singleton
class Method:
    identifier = "CombinedWig"
    menu_name = "A new Combined Wig file"
    usage_string = f"""
        Usage:
            {console_tools.subcommand_prefix} export combined_wig <comma-separated .wig files> <annotation_file> <output_file> [Optional Arguments]
            
        Optional Arguments:
            --n <string>     :=  Normalization method. Default: --n TTR
    """.replace("\n    ","\n")
    
    valid_cli_flags = [
        "--n",
        "--template",
        "-no-template",
    ]
    
    old_file_header = "#File: "
    
    @staticmethod
    @cli.add_command("export", "combined_wig")
    def from_args(args, kwargs):
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=3)
        
        Method.output(
            wig_list=args[0].split(","),
            annotation_path=args[1],
            output_path=args[2],
            normalization=kwargs["n"],
            metadata_path=kwargs["--template"] or (kwargs['no-template'] and ""), # None=default path, empty string=no template at all
        )
    
    @staticmethod
    def file_is_combined_wig(filepath):
        import os
        if os.path.exists(filepath):
            with open(filepath,'r') as f:
                for each in f.readlines():
                    if not each.startswith("#"):
                        return False
                    elif each.startswith("#"+Method.identifier):
                        return True
                    elif each.startswith(Method.old_file_header):
                        return True
                    elif each.startswith(Method.old_file_header.lower()):
                        return True
        return False
        
    @gui.add_menu("Pre-Processing", "Generate", menu_name)
    def on_menu_click(event):
        from pytransit.components import pop_up, panel_helpers
        
        @pop_up.create_pop_up(gui.frame, min_width=300)
        def create_pop_up_contents(pop_up_panel, sizer, refresh, close):
            normalization_getter   = panel_helpers.create_normalization_input(pop_up_panel, sizer)
            annotation_path_getter = panel_helpers.create_file_input(pop_up_panel, sizer, button_label="Select Annotation File", tooltip_text="", popup_title="Annotation File", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*', after_select=refresh)
            wig_paths_getter       = panel_helpers.create_multi_file_input(pop_up_panel, sizer, button_label="Select Wig Files", tooltip_text="", popup_title="Wig Files"      , default_folder=None, default_file_name="", allowed_extensions='Common output extensions (*.wig,*.tsv,*.dat,*.out)|*.wig;*.tsv;*.dat;*.out;|\nAll files (*.*)|*.*', after_select=refresh)
            
            @panel_helpers.create_button(pop_up_panel, sizer, label="Export")
            def when_button_clicked(event):
                normalization = normalization_getter()
                annotation_path = annotation_path_getter()
                wig_paths = wig_paths_getter()
                output_path = gui_tools.ask_for_output_file_path(
                    default_file_name=f"recent_export.comwig.tsv",
                    output_extensions='Common output extensions (*.comwig.tsv,*.tsv,*.csv,*.dat,*.out)|*.comwig.tsv;*.tsv;*.csv;*.dat;*.out;|\nAll files (*.*)|*.*',
                )
                
                # TODO: add validation here
                
                Method.output(
                    wig_list=wig_paths,
                    annotation_path=annotation_path,
                    output_path=output_path,
                    normalization=normalization,
                )
                
                logging.log(f"Finished ComWig Export: {output_path}")
                close()
    
    @staticmethod
    def output(*, annotation_path, wig_list, metadata_path=None, output_path=None, normalization=None, disable_logging=False):
        from pytransit.components.parameter_panel import progress_update
        # Defaults (even if argument directly provided as None)
        normalization     = normalization     if normalization     is not None else "TTR"
        output_path       = output_path       if output_path       is not None else "comwig.tsv"
        #metadata_path     = metadata_path     if metadata_path     is not None else os.path.dirname(output_path)+"/"+os.path.basename(output_path).split(".")[0]+".metadata.tsv"
        if metadata_path is None: metadata_path = "metadata.tsv" # justput in local dir, rather than trying to put in same dir as combined_wig file
        
        with transit_tools.TimerAndOutputs(method_name=Method.identifier, output_paths=[output_path], disable=disable_logging) as timer:
            logging.log("Starting Combined Wig Export")
            extra_info = LazyDict()
            
            # 
            # Read Data
            # 
            logging.log("Getting Data")
            read_counts_at_position_for_wig, ta_site_positions = tnseq_tools.CombinedWig.gather_wig_data(wig_list)
            ta_site_positions = ta_site_positions.astype(int)
            
            # 
            # normlize
            # 
            logging.log("Normalizing")
            (read_counts_at_position_for_wig, factors) = norm_tools.normalize_data(
                read_counts_at_position_for_wig,
                method=normalization,
                wig_list=wig_list,
                annotation_path=annotation_path,
            )
            extra_info.normalization_factors = misc.flatten_once(factors)
            
            # 
            # Genome Data
            # 
            positional_hash = transit_tools.get_pos_hash(annotation_path)
            rv2info = tnseq_tools.AnnotationFile(path=annotation_path).orf_to_info
            def get_gene_name(orf_id):
                return rv2info.get(orf_id,["-"])[0]
            
            # get ref genome from first wig file (we could check that they are all the same)
            # by this point, we know all the wig files exist and have same number of TA sites
            # this assumes 2nd line of each wig file is "variableStep chrom=XXX"; if not, set ref to "unknown"
            ref = "unknown"
            wig = wig_list[0]
            with open(wig) as file:
                line = file.readline()
                line = file.readline()
                if line.startswith("variableStep"):
                    splits = line.rstrip().split("=")
                    if len(splits)>=2: ref = splits[1]

            extra_info.reference_genome = ref
            number_of_ta_sites = len(ta_site_positions)
            
            # 
            # Compute rows
            # 
            rows = []
            for position_index, _ in enumerate(informative_iterator.ProgressBar(ta_site_positions, title="Running Export Method")):
                orf_info_sources = positional_hash.get(ta_site_positions[position_index], [])
                orf_ids          = [ orf for orf in orf_info_sources ]
                gene_names       = [ get_gene_name(orf) for orf in orf_info_sources ]
                read_counts      = [
                    # no decimals if raw counts
                    ("%d" if normalization == 'nonorm' else "%1.1f") % count 
                        for count in read_counts_at_position_for_wig[:,position_index]
                ]
                
                rows.append([ 
                    int(ta_site_positions[position_index]),
                    *read_counts,
                    *misc.flatten_once(zip(orf_ids, gene_names)),
                ])
                
                # Update progress
                if gui.is_active:
                    percentage = (100.0*position_index/number_of_ta_sites)
                    text = "Running Export Method... %5.1f%%" % percentage
                    if position_index % 1000 == 0:
                        progress_update(text, percentage)
            
            # 
            # wig fingerprints
            # 
            wig_fingerprints = wig_list
            
            #
            # Write main file
            #
            transit_tools.write_result(
                path=output_path,
                file_kind=Method.identifier,
                rows=rows,
                comments=[
                    f"Genome: {os.path.basename(annotation_path)}",
                    *[ f"File: {each}" for each in wig_fingerprints ],
                ],
                column_names=[
                    "TA Site Position",
                    *wig_fingerprints,
                    "ORF",
                    "Gene Name"
                ],
                extra_info={
                    **dict(
                        normalization=normalization,
                        annotation=os.path.basename(annotation_path),
                        wig_fingerprints=wig_fingerprints,
                    ),
                    **extra_info,
                },
            )
            logging.log("")  # Printing empty line to flush stdout
            logging.log("Finished Export")
            
            # 
            # metadata template
            # 
            if metadata_path:
                #wig_id_suggestions = tnseq_tools.wig_id_suggestions_from_filepaths(wig_list)
                wig_id_suggestions = [os.path.splitext(os.path.basename(x))[0].replace(' ','_') for x in wig_list] # strip off path and extension
                condition_name_suggestions = [ "CONDITION_NAME_HERE" ] * len(wig_id_suggestions)
                csv.write(
                    path=metadata_path,
                    column_names=["Id", "Condition", "Filename" ],
                    rows=list(zip(wig_id_suggestions, condition_name_suggestions, wig_fingerprints)),
                    seperator="\t",
                )
                logging.log(f"NOTE:\n    Edit the metadata file at: {metadata_path}")

    # 
    # Samples-area button
    # 
    @staticmethod
    @gui.add_wig_area_dropdown_option(name="Show Table")
    def click_show_table(event):
        from pytransit.components.spreadsheet import SpreadSheet
        selected_wigs = gui.selected_samples or gui.samples
        
        # 
        # heading (only if single wig)
        # 
        heading = ""
        if len(selected_wigs) == 1:
            heading = misc.human_readable_data(selected_wigs[0].extra_data)
        
        # 
        # row data
        # 
        column_names = [ 
            "TA Site Position",
            *[
                each.id 
                    for each in selected_wigs
            ] 
        ]
        ta_site_positions = selected_wigs[0].ta_sites
        rows = []
        for row_data in zip(*([ ta_site_positions ] + [ each.insertion_counts for each in selected_wigs ])):
            rows.append({
                column_name: cell_value
                    for column_name, cell_value in zip(column_names, row_data)
            })
        
        SpreadSheet(
            title="Read Counts",
            heading=heading,
            column_names=column_names,
            rows=rows,
            sort_by=[]
        ).Show()
    
    @staticmethod
    @gui.add_condition_area_dropdown_option(name="Show Table")
    def click_show_table(event):
        selected_wigs = gui.wigs_in_selected_conditions or gui.samples
        
        # 
        # row data
        # 
        column_names = [ 
            "TA Site Position",
            *[
                each.id 
                    for each in selected_wigs
            ] 
        ]
        ta_site_positions = selected_wigs[0].ta_sites
        rows = []
        for row_data in zip(*([ ta_site_positions ] + [ each.insertion_counts for each in selected_wigs ])):
            rows.append({
                column_name: cell_value
                    for column_name, cell_value in zip(column_names, row_data)
            })
        
        SpreadSheet(
            title="Read Counts",
            heading="",
            column_names=column_names,
            rows=rows,
            sort_by=[]
        ).Show()
