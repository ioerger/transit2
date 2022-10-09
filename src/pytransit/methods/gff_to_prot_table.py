from pytransit.components.parameter_panel import panel, progress_update
import sys
import os
import time

from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import misc, informative_iterator
from pytransit.globals import gui, cli, root_folder, debugging_enabled

@misc.singleton
class Method:
    name = "gff_to_prot"
    usage_string = f"""{console_tools.subcommand_prefix} convert gff_to_prot_table <annotation in gff format> <output file>"""
    
    inputs = LazyDict(
        annotation_path=None,
        output_file=None,
    )
    
    @cli.add_command("convert", name.lower())
    @classmethod
    def from_args(args, kwargs):
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=2)

        Method.inputs.update(dict(
            annotation_path= args[0],
            output_file=open(args[1], "w"),
        ))
        
        Method.Run()

    def Run(self):
        gff_file = open(self.inputs.annotation_path)
        output_file = self.inputs.output_file
        writer = csv.writer(output_file, delimiter="\t")
        lines = gff_file.readlines()
        gff_file.close()
        logging.log("Converting annotation file from GFF3 format to prot_table format")

        for i, line in enumerate(lines):
            line = line.strip()
            if len(line) == 0 or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                sys.stderr.write(("Ignoring invalid row with entries: {0}\n".format(cols)))
                continue
            if (cols[2]) == "CDS":  # if you also want tRNAs and rRNAs, modify here
                if "locus_tag" not in line:
                    print("warning: skipping lines that do not contain 'locus_tag'")
                    continue
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6].strip()
                size = int(abs(end - start + 1) / 3)  # includes stop codon
                labels = {}
                for pair in cols[8].split(";"):
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
        for pair in cols[8].split(";"):
            k, v = pair.split("=")
            labels[k] = v

        if (cols[2]) == "CDS" and labels["Parent"] == parent:
            # return labels.get("Note", '-')
            return labels.get("product", "-")
        return "-"