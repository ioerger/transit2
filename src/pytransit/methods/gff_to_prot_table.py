from pytransit.components.parameter_panel import panel, progress_update
import sys
import os
import time

from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.specific_tools.tnseq_tools import GffFile
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import misc, informative_iterator
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
        
@misc.singleton
class Method:
    name = "gff_to_prot"
    menu_name = "Prot Table from GFF file"
    usage_string = f"""Usage: {console_tools.subcommand_prefix} convert gff_to_prot_table <gff_file> <output_file>"""
    
    inputs = LazyDict(
        path_to_gff=None,
        output_file=None,
    )
    
    @staticmethod
    @cli.add_command("convert", name.lower())
    def from_args(args, kwargs):
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=2)

        Method.inputs.update(dict(
            path_to_gff= args[0],
            output_file=open(args[1], "w"),
        ))
        
        Method.Run()
    
    @gui.add_menu("Pre-Processing", "Generate", menu_name)
    def on_menu_click(event):
        from pytransit.components import pop_up
        from pytransit.components import panel_helpers
        
        @pop_up.create_pop_up(gui.frame, min_width=300)
        def create_pop_up_contents(pop_up_panel, sizer, refresh, close):
            gff_path_getter = panel_helpers.create_file_input(pop_up_panel, sizer, button_label="Add GFF File", tooltip_text="", popup_title="GFF File", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*', after_select=refresh)
            
            @panel_helpers.create_button(pop_up_panel, sizer, label="Convert")
            def when_button_clicked(event):
                import os
                path_to_gff = gff_path_getter()
                # get name
                name = os.path.basename(path_to_gff).replace(".gff3","").replace(".gff","")
                # get output path
                output_path = gui_tools.ask_for_output_file_path(
                    default_file_name=f"{name}.prot_table",
                    output_extensions=['prot_table', 'tsv', 'csv', 'dat', 'txt'],
                )
                Method.inputs.update(dict(
                    path_to_gff=path_to_gff,
                    output_file=open(output_path, "w"),
                ))
                Method.Run()
                logging.log(f"Conversion complete, written to {output_path}")
    
    def Run(self):
        gff = GffFile(path=self.inputs.path_to_gff)
        import csv
        output_file = self.inputs.output_file
        writer = csv.writer(output_file, delimiter="\t")
        logging.log("Converting annotation file from GFF3 format to prot_table format")
        
        for row in gff.as_prot_table_rows():
            writer.writerow(row)
        output_file.close()
        logging.log("Finished conversion")