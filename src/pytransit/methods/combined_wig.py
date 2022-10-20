from pytransit.components.parameter_panel import panel, progress_update
import sys
import os
import time

import numpy

from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import misc, informative_iterator
from pytransit.components import samples_area
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    identifier = "CombinedWig"
    menu_name = "Combined Wig"
    usage_string = f"""
        {console_tools.subcommand_prefix} export combined_wig <comma-separated .wig files> <annotation .prot_table> <output file> [-n normalization_method]
        
            default normalization_method=TTR
    """.replace("\n    ","\n")
    
    inputs = LazyDict(
        ctrldata=None,
        normalization=None,
        annotation_path=None,
        output_path="latest.comwig.csv",
        ref=None,
    )
    
    @cli.add_command("export", "combined_wig")
    @staticmethod
    def from_args(args, kwargs):
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=3)

        Method.inputs.update(dict(
            ctrldata= args[0].split(","),
            annotation_path=args[1],
            output_path=args[2],
            normalization= kwargs.get("n", "TTR"),
        ))
        
        Method.Run()
        
    @gui.add_menu("Pre-Processing", "Export", menu_name)
    def on_menu_click(event):
        from pytransit.components import pop_up, panel_helpers
        
        @pop_up.create_pop_up(gui.frame, min_width=300)
        def create_pop_up_contents(pop_up_panel, sizer, refresh, close):
            normalization_getter   = panel_helpers.create_normalization_input(pop_up_panel, sizer)
            annotation_path_getter = panel_helpers.create_file_input(pop_up_panel, sizer, button_label="Select Annotation File", tooltip_text="", popup_title="Annotation File", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*', after_select=refresh)
            wig_paths_getter       = panel_helpers.create_multi_file_input(pop_up_panel, sizer, button_label="Select Wig Files", tooltip_text="", popup_title="Wig Files"      , default_folder=None, default_file_name="", allowed_extensions='Common output extensions (*.wig,*.csv,*.dat,*.out)|*.wig;*.csv;*.dat;*.out;|\nAll files (*.*)|*.*', after_select=refresh)
            
            @panel_helpers.create_button(pop_up_panel, sizer, label="Export")
            def when_button_clicked(event):
                normalization = normalization_getter()
                annotation_path = annotation_path_getter()
                wig_paths = wig_paths_getter()
                output_path = gui_tools.ask_for_output_file_path(
                    default_file_name=f"recent_export.comwig.csv",
                    output_extensions='Common output extensions (*.comwig.csv,*.csv,*.dat,*.out)|*.comwig.csv;*.csv;*.dat;*.out;|\nAll files (*.*)|*.*',
                )
                print(f'''wig_paths = {wig_paths}''')
                
                # TODO: add validation here
                
                Method.inputs.update(dict(
                    ctrldata=wig_paths,
                    annotation_path=annotation_path,
                    output_path=output_path,
                    normalization=normalization,
                ))
                
                Method.Run()
                logging.log(f"Finished ComWig Export: {output_path}")
                close()
       
    
    def Run(self):
        logging.log("Starting Combined Wig Export")
        start_time = time.time()
        extra_info = LazyDict()
        
        # 
        # Read Data
        # 
        logging.log("Getting Data")
        read_counts_at_position_for_wig, ta_site_positions = tnseq_tools.CombinedWig.gather_wig_data(self.inputs.ctrldata)
        ta_site_positions = ta_site_positions.astype(int)
        
        # 
        # normlize
        # 
        logging.log("Normalizing")
        (read_counts_at_position_for_wig, factors) = norm_tools.normalize_data(
            read_counts_at_position_for_wig,
            method=self.inputs.normalization,
            wig_list=self.inputs.ctrldata,
            annotation_path=self.inputs.annotation_path,
        )
        if type(factors[0]) == type(0.0):
            extra_info.normalization_factors = factors.flatten().tolist()
        else:
            extra_info.normalization_factors = misc.flatten_once(factors.tolist())
        
        # 
        # Genome Data
        # 
        positional_hash = transit_tools.get_pos_hash(self.inputs.annotation_path)
        rv2info = transit_tools.get_gene_info(self.inputs.annotation_path)
        def get_gene_name(orf_id):
            return rv2info.get(orf_id,["-"])[0]
        
        # get ref genome from first wig file (we could check that they are all the same)
        # by this point, we know all the wig files exist and have same number of TA sites
        # this assumes 2nd line of each wig file is "variableStep chrom=XXX"; if not, set ref to "unknown"
        self.inputs.ref = "unknown"
        wig = self.inputs.ctrldata[0]
        with open(wig) as file:
            line = file.readline()
            line = file.readline()
            if line.startswith("variableStep"):
                splits = line.rstrip().split("=")
                if len(splits)>=2: self.inputs.ref = splits[1]

        extra_info.reference_genome = self.inputs.ref
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
                ("%d" if self.inputs.normalization == 'nonorm' else "%1.1f") % count 
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
        wig_fingerprints = tnseq_tools.filepaths_to_fingerprints(self.inputs.ctrldata)
        
        # 
        # Stats per wig
        # 
        from pytransit.methods.tnseq_stats import Method as tnseq_stats
        logging.log("Generating stats")
        stats = LazyDict()
        for wig_file_index, each_wig_fingerprint in enumerate(wig_fingerprints):
            wig_insertion_counts = read_counts_at_position_for_wig[wig_file_index, :]
            
            (
                density,
                mean_read,
                non_zero_mean_read,
                non_zero_median_read,
                max_read,
                total_read,
                skew,
                kurtosis,
            ) = tnseq_tools.get_data_stats(wig_insertion_counts)
            
            import ez_yaml
            import sys
            ez_yaml.yaml.version = None
            ez_yaml.yaml.width = sys.maxint if hasattr(sys, "maxint") else sys.maxsize
            stats[each_wig_fingerprint] = ez_yaml.to_string(dict(
                density=float(density),
                mean_read=float(mean_read),
                non_zero_mean_read=float(non_zero_mean_read),
                non_zero_median_read= 0 if numpy.isnan(non_zero_median_read) else int(non_zero_median_read),
                max_read=int(max_read),
                total_read=int(total_read),
                skew=float(skew),
                kurtosis=float(kurtosis),
                pickands_tail_index=float(tnseq_stats.pickands_tail_index(wig_insertion_counts)),
            )).replace("\n", ", ").strip()
        
        #
        # Write to output
        #
        transit_tools.write_result(
            path=self.inputs.output_path,
            file_kind=Method.identifier,
            rows=rows,
            column_names=[
                "TA Site Position",
                *wig_fingerprints,
                "ORF",
                "Gene Name"
            ],
            extra_info={
                **dict(
                    time=(time.time() - start_time),
                    normalization=self.inputs.normalization,
                    annotation=os.path.basename(self.inputs.annotation_path),
                    wig_fingerprints=wig_fingerprints,
                ),
                **extra_info,
                "stats": dict(stats),
            },
        )

        logging.log("")  # Printing empty line to flush stdout
        logging.log("Finished Export")
        
    # 
    # Samples-area button
    # 
    @gui.add_wig_area_dropdown_option(name="Show Table")
    @staticmethod
    def click_show_table(event):
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
        ta_site_positions = selected_wigs[0].positions
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
    
    @gui.add_condition_area_dropdown_option(name="Show Table")
    @staticmethod
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
        ta_site_positions = selected_wigs[0].positions
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