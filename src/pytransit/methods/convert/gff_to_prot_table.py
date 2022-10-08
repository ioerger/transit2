import sys
import os

from pytransit.tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub

import csv
import traceback
from pytransit.tools import transit_tools, console_tools, logging
from pytransit.methods import convert_base as base
from pytransit.globals import gui, cli, root_folder, debugging_enabled

############# Description ##################

short_name = "gff_to_prot"
long_name = "GFF to Prot_table converter"
description = "Convert a GFF file to Prot_table format"
label = "GFF3 to Prot_table"

############# Analysis Method ##############


class Converter(base.TransitConvert):
    def __init__(self):
        base.TransitConvert.__init__(
            self, short_name, long_name, description, label, GffProtMethod, GffProtGUI
        )


################# GUI ##################
class GffProtGUI(base.ConvertGUI):
    def __init__(self):
        base.ConvertGUI.__init__(self)


########## METHOD #######################


class GffProtMethod(base.ConvertMethod):
    """
    GffProtMethod
    """

    def __init__(self, annotation_path, output, wxobj=None):
        self.short_name = short_name
        self.long_name = long_name
        self.description = description
        self.label = label
        self.output = output
        self.annotation_path = annotation_path
        base.ConvertMethod.__init__(
            self,
            short_name,
            long_name,
            description,
            label,
            annotation_path,
            output,
            wxobj=wxobj,
        )

    @classmethod
    def from_gui(cls, wxobj):
        """ """
        # Get Annotation file
        from pytransit.globals import gui, cli, root_folder, debugging_enabled
        annotation_path = gui.annotation_path
        if not transit_tools.validate_annotation(annotation_path):
            return None

        # Get output path
        defaultFileName = "{0}.prot_table".format(
            os.path.splitext(os.path.basename(annotation_path))[0]
        )
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path:
            return None
        output_file = open(output_path, "w")

        return cls(annotation_path, output_file, wxobj)
    
    @cli.add_command("convert", "gff_to_prot")
    @classmethod
    def from_args(self, args, kwargs):
        console_tools.enforce_number_of_args(args, self.usage_string, at_least=2)

        annotation_path = args[0]
        outpath = args[1]
        output_file = open(outpath, "w")

        return self(annotation_path, output_file)

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

    def Run(self):
        gff_file = open(self.annotation_path)
        output_file = self.output
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

    usage_string = f"""{console_tools.subcommand_prefix} convert gff_to_prot_table <annotation in gff format> <output file>"""
