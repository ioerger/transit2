from pytransit.components.parameter_panel import panel, progress_update
import sys
import os
import time

from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.specific_tools.tnseq_tools import GffFile
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import misc, informative_iterator
from pytransit.globals import gui, cli, root_folder, debugging_enabled
        
@misc.singleton
class Method:
    name = "gff_to_prot"
    menu_name = "GFF to Prot Table"
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
    
    @gui.add_menu("Pre-Processing", "Convert", menu_name)
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
        gff_file = open(self.inputs.path_to_gff)
        output_file = self.inputs.output_file
        import csv
        writer = csv.writer(output_file, delimiter="\t")
        lines = gff_file.readlines()
        gff_file.close()
        logging.log("Converting annotation file from GFF3 format to prot_table format")
        
        # FIXME: should probably use transit_tools.GffFile to load/parse the data
        for i, line in enumerate(lines):
            line = line.strip()
            if len(line) == 0 or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < GffFile.max_number_of_columns:
                sys.stderr.write(("Ignoring invalid row with entries: {0}\n".format(cols)))
                continue
            if (cols[GffFile.end_index]) == "CDS":  # if you also want tRNAs and rRNAs, modify here
                if "locus_tag" not in line:
                    print("warning: skipping lines that do not contain 'locus_tag'")
                    continue
                start = int(cols[GffFile.start_index])
                end = int(cols[GffFile.end_index])
                strand = cols[GffFile.strand_index].strip()
                size = int(abs(end - start + 1) / 3)  # FIXME: why divide by 3?
                labels = {}
                for pair in cols[GffFile.attributes_index].split(";"):
                    k, v = pair.split("=")
                    labels[k.strip()] = v.strip()
                Rv = labels["locus_tag"].strip()  # error out if not found
                gene = labels.get("gene", "")  # or Name?
                if gene == "":
                    gene = "-"
                desc = labels.get("product", "")
                vals = [desc, start, end, strand, size, "-", "-", gene, Rv, "-"]
                writer.writerow(vals)
        output_file.close()
        logging.log("Finished conversion")
    
    def get_description(self, line, parent):
        cols = line.split("\t")
        labels = {}
        print(line)
        print(len(cols))
        for pair in cols[GffFile.attributes_index].split(";"):
            k, v = pair.split("=")
            labels[k] = v

        if (cols[GffFile.end_index]) == "CDS" and labels["Parent"] == parent:
            # return labels.get("Note", '-')
            return labels.get("product", "-")
        return "-"


# UNUSED at the moment
def annotation_gff3_to_pt(event):
    with gui_tools.nice_error_log:
        path_to_gff = gui.path_to_gff
        default_file = transit_tools.fetch_name(path_to_gff) + ".prot_table"
        # default_dir = os.path.dirname(os.path.realpath(__file__))
        default_dir = os.getcwd()

        if not path_to_gff:
            # NOTE: was a popup
            logging.error("Error: No annotation file selected.")
        else:
            output_path = frame.SaveFile(default_dir, default_file)
            if not output_path:
                return
            logging.log("Converting annotation file from GFF3 format to prot_table format")

            output = open(output_path, "w")
            with open(path_to_gff) as file:
                for line in file:
                    if line.startswith("#"):
                        continue
                    tmp = line.strip().split("\t")
                    chr = tmp[0]
                    type = tmp[GffFile.type_index]
                    start = int(tmp[GffFile.start_index])
                    end = int(tmp[GffFile.end_index])
                    length = ((end - start + 1) / GffFile.start_index) - 1
                    strand = tmp[GffFile.strand_index]
                    features = dict([tuple(f.split("=")) for f in tmp[GffFile.attributes_index].split(";")])
                    if "ID" not in features:
                        continue
                    orf = features["ID"]
                    name = features.get("Name", features.get("Gene Name","-"))
                    if name == "-":
                        name = features.get("name", "-")

                    desc = features.get("Description", "-")
                    if desc == "-":
                        desc = features.get("description", "-")
                    if desc == "-":
                        desc = features.get("Desc", "-")
                    if desc == "-":
                        desc = features.get("desc", "-")

                    someID = "-"
                    someID2 = "-"
                    COG = "-"
                    output.write(
                        "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n"
                        % (
                            desc,
                            start,
                            end,
                            strand,
                            length,
                            someID,
                            someID2,
                            name,
                            orf,
                            COG,
                        )
                    )
            output.close()
            logging.log("Finished conversion")