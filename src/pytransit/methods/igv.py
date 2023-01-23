from pytransit.components.parameter_panel import panel, progress_update
import sys
import os
import time

from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import misc, informative_iterator
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

@misc.singleton
class Method:
    name = "igv"
    description = "A method to export datasets in 'IGV' format."
    menu_name = f"IGV file from wig file"
    usage_string = f"""
        Usage:
            {console_tools.subcommand_prefix} export igv <comma-separated .wig files> <annotation_file> <output_file>
    """.replace("\n        ", "\n")
    
    inputs = LazyDict(
        ctrldata=None,
        normalization=None,
        annotation_path=None,
        output_file=None,
    )
    
    @staticmethod
    @cli.add_command("export", name.lower())
    def from_args(args, kwargs):
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=3)

        Method.inputs.update(dict(
            ctrldata= args[0].split(","),
            annotation_path= args[1],
            output_file=open(args[2], "w"),
            normalization= kwargs.get("n", "TTR"),
        ))
        
        Method.Run()

    @gui.add_menu("Pre-Processing", "Generate", menu_name)
    def on_menu_click(event):
        from pytransit.components import pop_up
        from pytransit.components import panel_helpers
        
        @pop_up.create_pop_up(gui.frame, min_width=300)
        def create_pop_up_contents(pop_up_panel, sizer, refresh, close):
            # TODO: consider allowing selecting a wig-id from a dropdown
            # wig_obj_getter         = panel_helpers.create_wig_choice(pop_up_panel, sizer, label_text="Select Wig Id")
            normalization_getter   = panel_helpers.create_normalization_input(pop_up_panel, sizer)
            wig_path_getter        = panel_helpers.create_file_input(pop_up_panel, sizer, button_label="Select Wig File",        tooltip_text="", popup_title="Wig File"       , default_folder=None, default_file_name="", allowed_extensions=['wig','tsv','csv','dat','out'], after_select=refresh)
            annotation_path_getter = panel_helpers.create_file_input(pop_up_panel, sizer, button_label="Select Annotation File", tooltip_text="", popup_title="Annotation File", default_folder=None, default_file_name="", after_select=refresh)
            
            @panel_helpers.create_button(pop_up_panel, sizer, label="Export")
            def when_button_clicked(event):
                wig_path = wig_path_getter()
                normalization = normalization_getter()
                annotation_path = annotation_path_getter()
                output_path = gui_tools.ask_for_output_file_path(
                    default_file_name=f"recent_export.igv",
                    output_extensions=['igv', 'tsv', 'csv', 'dat', 'txt'],
                )
                
                # TODO: add validation here
                
                Method.inputs.update(dict(
                    ctrldata=wig_path,
                    annotation_path=annotation_path,
                    output_file=open(output_path, "w"),
                    normalization=normalization,
                ))
                
                Method.Run()
                logging.log(f"Finished IGV Export: {output_path}")
    
    def Run(self):
        logging.log("Starting IGV Export")
        start_time = time.time()

        # Get orf data
        logging.log("Getting Data")
        (fulldata, position) = tnseq_tools.CombinedWig.gather_wig_data(self.inputs.ctrldata)
        (fulldata, factors) = norm_tools.normalize_data(
            fulldata, self.inputs.normalization, self.inputs.ctrldata, self.inputs.annotation_path
        )
        position = position.astype(int)

        hash = transit_tools.get_pos_hash(self.inputs.annotation_path)
        rv2info = tnseq_tools.AnnotationFile(path=self.inputs.annotation_path).orf_to_info

        logging.log("Normalizing")
        self.inputs.output_file.write("#Converted to IGV with TRANSIT.\n")
        if self.inputs.normalization != "nonorm":
            self.inputs.output_file.write("#Reads normalized using '%s'\n" % self.inputs.normalization)
            if type(factors[0]) == type(0.0):
                self.inputs.output_file.write(
                    "#Normalization Factors: %s\n"
                    % "\t".join(["%s" % f for f in factors.flatten()])
                )
            else:
                self.inputs.output_file.write(
                    "#Normalization Factors: %s\n"
                    % " ".join([",".join(["%s" % bx for bx in b]) for b in factors])
                )

        self.inputs.output_file.write("#Files:\n")
        for f in self.inputs.ctrldata:
            self.inputs.output_file.write("#%s\n" % f)

        dataset_str = "\t".join([transit_tools.fetch_name(F) for F in self.inputs.ctrldata])
        self.inputs.output_file.write("#Chromosome\tStart\tEnd\tFeature\t%s\tTAs\n" % dataset_str)
        chrom = transit_tools.fetch_name(self.inputs.annotation_path)

        (K, N) = fulldata.shape
        
        for i, pos in enumerate(position):
            self.inputs.output_file.write(
                "%s\t%s\t%s\tTA%s\t%s\t1\n"
                % (
                    chrom,
                    position[i],
                    position[i] + 1,
                    position[i],
                    "\t".join(["%1.1f" % fulldata[j][i] for j in range(len(fulldata))]),
                )
            )

            # Update progress
            percentage = (100.0 * i / N)
            text = "Running Export Method... %5.1f%%" % percentage
            progress_update(text, percentage)
        self.inputs.output_file.close()

        logging.log("")  # Printing empty line to flush stdout
        logging.log("Finished Export")