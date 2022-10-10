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
    name = "mean_counts"
    description = "A method to export and normalized datasets in 'Mean Gene Counts' format."
    usage_string = f"""{console_tools.subcommand_prefix} export mean_counts <comma-separated .wig files>|<combined_wig> <annotation .prot_table> <output file> [-c]\n note: append -c if inputing a combined_wig file\n"""
    
    inputs = LazyDict(
        combined_wig=None,
        ctrldata=None,
        normalization=None,
        annotation_path=None,
        output_file=None,
    )
    
    @cli.add_command("export", name.lower())
    @staticmethod
    def from_args(args, kwargs):
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=3)

        Method.inputs.update(dict(
            combined_wig=kwargs.get("c", False),
            ctrldata=args[0].split(","),
            annotation_path=args[1],
            output_file=open(args[2], "w"),
            normalization=kwargs.get("n", "TTR"),
        ))
        
        Method.Run()

    def Run(self):
        logging.log("Starting Gene Mean Counts Export")
        start_time = time.time()

        # Get orf data
        logging.log("Getting Data")
        if self.inputs.combined_wig:
            (position, fulldata, datasets) = tnseq_tools.CombinedWigData.load(
                self.inputs.ctrldata[0]
            )
        else:
            (fulldata, position) = tnseq_tools.CombinedWig.gather_wig_data(self.inputs.ctrldata)
        (fulldata, factors) = norm_tools.normalize_data(
            fulldata, self.inputs.normalization, self.inputs.ctrldata, self.inputs.annotation_path
        )
        position = position.astype(int)

        hash = transit_tools.get_pos_hash(self.inputs.annotation_path)
        rv2info = transit_tools.get_gene_info(self.inputs.annotation_path)

        logging.log("Normalizing")
        self.inputs.output_file.write("#Summarized to Mean Gene Counts with TRANSIT.\n")
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
        names = datasets if self.inputs.combined_wig else self.inputs.ctrldata
        for f in names:
            self.inputs.output_file.write("#%s\n" % f)

        K, Nsites = fulldata.shape
        # Get Gene objects
        if self.inputs.combined_wig:
            G = tnseq_tools.Genes(
                self.inputs.ctrldata,
                self.inputs.annotation_path,
                norm=self.inputs.normalization,
                data=fulldata,
                position=position,
            )
        else:
            G = tnseq_tools.Genes(
                self.inputs.ctrldata, self.inputs.annotation_path, norm=self.inputs.normalization
            )
        N = len(G)
        
        if self.inputs.combined_wig:
            dataset_header = "\t".join(datasets)
        else:
            dataset_header = "\t".join(
                [transit_tools.fetch_name(D) for D in self.inputs.ctrldata]
            )
        self.inputs.output_file.write("#Orf\tName\tNumber of TA sites\t%s\n" % dataset_header)
        for i, gene in enumerate(G):
            if gene.n > 0:
                data_str = "\t".join(["%1.2f" % (M) for M in numpy.mean(gene.reads, 1)])
            else:
                data_str = "\t".join(["%1.2f" % (Z) for Z in numpy.zeros(K)])
            self.inputs.output_file.write(
                "%s\t%s\t%s\t%s\n" % (gene.orf, gene.name, gene.n, data_str)
            )

            # Update progress
            percentage = (100.0 * i / N)
            text = "Running Export Method... %5.1f%%" % percentage
            progress_update(text, percentage)
        self.inputs.output_file.close()

        logging.log("")  # Printing empty line to flush stdout
        logging.log("Finished Export")
