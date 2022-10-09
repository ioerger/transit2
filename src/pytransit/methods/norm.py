from pytransit.components.parameter_panel import panel, progress_update
import sys
import os
import time

from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools, informative_iterator
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import misc
from pytransit.globals import gui, cli, root_folder, debugging_enabled

@misc.singleton
class Method:
    name = "normalize"
    description = "Method for normalizing datasets and outputting into CombinedWig file."
    usage_string = f"""
        {console_tools.subcommand_prefix} norm <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    
            Optional Arguments:
            -n <string>     :=  Normalization method. Default: -n TTR
    """.replace("\n        ","\n")
    
    
    # 
    # newer method
    # 
    @cli.add_command("normalize")
    @staticmethod
    def from_args(args, kwargs):
        is_combined_wig = "c" in kwargs
        
        if not is_combined_wig:
            console_tools.enforce_number_of_args(args, Method.usage_string, at_least=2)
            infile_path = kwargs.get("c")  # only 1 input wig file
            otuput_path = args[0]  # if no arg give, could print to screen
        else:
            console_tools.enforce_number_of_args(args, Method.usage_string, at_least=1)
            infile_path = args[0]  # only 1 input wig file
            otuput_path = args[1]  # if no arg give, could print to screen
        
        normalization = kwargs.get("n", "TTR")
        combined_wig = is_combined_wig
        
        Method.run_normalize(
            combined_wig=is_combined_wig,
            infile_path=infile_path,
            output_path=output_path,
            normalization=
        )

    @staticmethod
    def run_normalize(combined_wig, infile_path, output_path, normalization):
        logging.log("Starting Normalization")
        start_time = time.time()

        # determine ref genome from first; assume they are all the same; assume wigs have 2 header lines
        line2 = "variableStep chrom="  # unknown
        with open(infile_path) as file:
            for line in file:
                if line.startswith("variableStep"):
                    line2 = line.rstrip()
                    break

        if combined_wig == True:
            (sites, data, files) = tnseq_tools.read_combined_wig(infile_path)
        else:
            (data, sites) = tnseq_tools.CombinedWig.gather_wig_data(infile_path)
        (data, factors) = norm_tools.normalize_data(data, normalization)

        print("writing", output_path)
        with open(filepath,'w') as file:
            file.write("# %s normalization of %s\n" % (normalization, infile_path))
            if combined_wig == True:
                for each_wig_file in files:
                    file.write("#File: %s\n" % each_wig_file)
                for i in range(len(sites)):
                    file.write(
                        "\t".join(
                            [str(sites[i])] + ["%0.1f" % x for x in list(data[..., i])]
                        )
                        + "\n"
                    )
            else:
                file.write(line2 + "\n")
                for site_index in range(len(sites)):
                    file.write("%s %s\n" % (sites[site_index], int(data[0, site_index])))
        
        logging.log("Finished Normalization")

    
    # 
    # older method
    # 
    @cli.add_command("export", "norm")
    @staticmethod
    def from_export(args, kwargs):
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=3)
        Method.run_norm(
            ctrldata=args[0].split(","),
            annotation_path=args[1],
            output_path=args[2],
            normalization=kwargs.get("n", "TTR"),
        )
    
    @staticmethod
    def run_norm(ctrldata, annotation_path, output_path, normalization):
        logging.log("Starting Normalization")
        transit_tools.convert_to_combined_wig(
            ctrldata,
            annotation_path,
            output_path,
            normchoice=normalization,
        )
        logging.log("Finished Normalization")