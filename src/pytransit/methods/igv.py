from pytransit.components.parameter_panel import panel, progress_update
import sys
import os
import time

from pytransit.tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools, informative_iterator
from pytransit.basics.lazy_dict import LazyDict
from pytransit.basics import misc
from pytransit.globals import gui, cli, root_folder, debugging_enabled

@misc.singleton
class Method:
    name = "igv"
    description = "A method to export and normalized datasets in 'IGV' format."
    usage_string = f"""{console_tools.subcommand_prefix} export igv <comma-separated .wig files> <annotation .prot_table> <output file>"""
    
    inputs = LazyDict(
        ctrldata=None,
        normalization=None,
        annotation_path=None,
        output_file=None,
    )
    
    @cli.add_command("export", name.lower())
    @staticmethod
    def from_args(args, kwargs):
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=3)

        Method.inputs.update(dict(
            ctrldata= args[0].split(","),
            annotation_path= args[1],
            output_file=open(args[2], "w"),
            normalization= kwargs.get("n", "TTR"),
        ))
        
        Method.Run()

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
        rv2info = transit_tools.get_gene_info(self.inputs.annotation_path)

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