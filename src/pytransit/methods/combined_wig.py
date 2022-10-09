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
    usage_string = f"""
        {console_tools.subcommand_prefix} export combined_wig <comma-separated .wig files> <annotation .prot_table> <output file> [-n normalization_method]
        
            default normalization_method=TTR
    """.replace("\n    ","\n")
    
    inputs = LazyDict(
        ctrldata=None,
        normalization=None,
        annotation_path=None,
        output_file=None,
        ref=None,
    )
    
    @cli.add_command("export", "combined_wig")
    @staticmethod
    def from_args(args, kwargs):
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=3)

        Method.inputs.update(dict(
            ctrldata= args[0].split(","),
            annotation_path= args[1],
            output_file=open(args[2], "w"),
            normalization= kwargs.get("n", "TTR"),
        ))
        
        Method.Run()

    def Run(self):
        logging.log("Starting Combined Wig Export")
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
        self.inputs.output_file.write("#Converted to CombinedWig with TRANSIT.\n")
        self.inputs.output_file.write("#normalization method: %s\n" % self.inputs.normalization)
        if self.inputs.normalization != "nonorm":
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

        # get ref genome from first wig file (we could check that they are all the same)
        # by this point, we know all the wig files exist and have same number of TA sites
        # this assumes 2nd line of each wig file is "variableStep chrom=XXX"; if not, set ref to "unknown"
        self.inputs.ref = "unknown"
        wig = self.inputs.ctrldata[0]
        with open(wig) as file:
          line = file.readline()
          line = file.readline()
          if line.startswith("variableStep"): 
            w = line.rstrip().split("=")
            if len(w)>=2: self.inputs.ref = w[1]

        (K,N) = fulldata.shape
        self.inputs.output_file.write("#RefGenome: %s\n" % self.inputs.ref)
        for f in self.inputs.ctrldata:
            self.inputs.output_file.write("#File: %s\n" % f)
        self.inputs.output_file.write("#TA_coord\t%s\n" % ('\t'.join(self.inputs.ctrldata)))

        for i,pos in enumerate(position):
            #self.inputs.output_file.write("%d\t%s\t%s\n" % (position[i], "\t".join(["%1.1f" % c for c in fulldata[:,i]]),",".join(["%s (%s)" % (orf,rv2info.get(orf,["-"])[0]) for orf in hash.get(position[i], [])])   ))
            if self.inputs.normalization != 'nonorm':
                vals = "\t".join(["%1.1f" % c for c in fulldata[:,i]])
            else:
                # no decimals if raw counts
                vals = "\t".join(["%d"    % c for c in fulldata[:,i]])
            self.inputs.output_file.write(
                "%d\t%s\t%s\n" % (
                    position[i],
                    vals,
                    ",".join(
                        [
                            "%s (%s)" % (
                                orf,
                                rv2info.get(orf,["-"])[0]
                            )
                                for orf in hash.get(position[i], [])
                        ]
                    )
                )
            )
            # Update progress
            percentage = (100.0*i/N)
            text = "Running Export Method... %5.1f%%" % percentage
            if i % 1000 == 0:
                progress_update(text, percentage)
        self.inputs.output_file.close()

        logging.log("")  # Printing empty line to flush stdout
        logging.log("Finished Export")