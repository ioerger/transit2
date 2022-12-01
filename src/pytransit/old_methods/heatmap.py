import sys
import os
import time
import math
import random
import datetime

import numpy

from pytransit.old_methods import analysis_base as base
from pytransit.specific_tools.transit_tools import HAS_R, r, DataFrame, globalenv, IntVector, FloatVector, StrVector, rpackages
form pytransit.specific_tools import transit_tools, tnseq_tools, norm_tools, stat_tools, console_tools, logging


############# Description ##################

short_name = "heatmap"
long_name = "Heatmap"
short_desc = "Heatmap among Conditions"
long_desc = "Heatmap among Conditions"
transposons = ["himar1", "tn5"]

columns = ["Position", "Reads", "Genes"]  # ???

############# Analysis Method ##############


class Analysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            transposons,
            HeatmapMethod,
            HeatmapGUI,
            [],
        )


################## FILE ###################

# there is no output file that could be loaded into the GUI

# class HeatmapFile(base.TransitFile):
#
#    def __init__(self):
#        base.TransitFile.__init__(self, "#CombinedWig", columns)
#
#    def get_header(self, path):
#        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
#        return text

################# GUI ##################

# right now, tnseq_stats is just intended for the command-line; TRI


class HeatmapGUI(base.AnalysisGUI):
    def __init__(self):
        base.AnalysisGUI.__init__(self)


########## METHOD #######################

# should Heatmap be a SingleConditionMethod? args like normalization are irrelevant


class HeatmapMethod(base.SingleConditionMethod):
    usage_string = "usage: python3 %s heatmap <anova_or_zinb_output> <heatmap.png> -anova|-zinb [-topk <int>] [-qval <float>] [-low_mean_filter <int>]\n note: genes are selected based on qval<0.05 by default" % sys.argv[0]

    def __init__(self, gene_means, outfile):
        ctrldata = None  # initializers for superclass
        annotation_path = ""
        output_file = outfile
        replicates = "Sum"
        normalization = "nonorm"
        LOESS = False
        ignore_codon = True
        n_terminus = 0.0
        c_terminus = 0.0
        wxobj = None
        base.SingleConditionMethod.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            ctrldata,
            annotation_path,
            output_file,
            replicates=replicates,
            normalization=normalization,
            LOESS=LOESS,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            wxobj=wxobj,
        )

    @classmethod
    def from_args(self, args, kwargs):
        console_tools.enforce_number_of_args(args, self.usage_string, at_least=3)

        self.inputs.filetype = None
        if kwargs.get("anova", False):
            self.inputs.filetype = "anova"
        elif kwargs.get("zinb", False):
            self.inputs.filetype = "zinb"
        else:
            logging.error(f"requires --anova or --zinb argument, see usage string below.\n{self.usage_string}")
        
        self.inputs.input_path = args[0]
        self.inputs.output_path = args[1]
        self.inputs.adj_p_value = float(kwargs.get("qval", 0.05))
        self.inputs.top_k = int(kwargs.get("topk", -1))
        self.inputs.low_mean_filter = int(
            kwargs.get("low_mean_filter", 5)
        )  # filter out genes with grandmean<5 by default
        return self(self.inputs.input_path, outfile=self.inputs.output_path)

    def Run(self):
        transit_tools.require_r_to_be_installed()
        if self.inputs.filetype != "anova" and self.inputs.filetype != "zinb":
            logging.error("filetype not recognized: %s" % self.inputs.filetype)

        headers = None
        data, hits = [], []
        n = -1  # number of conditions

        with open(self.inputs.input_path) as file:
            for line in file:
                w = line.rstrip().split("\t")
                if line[0] == "#" or (
                    "pval" in line and "padj" in line
                ):  # check for 'pval' for backwards compatibility
                    headers = w
                    continue  # keep last comment line as headers
                # assume first non-comment line is header
                if n == -1:
                    # ANOVA header line has names of conditions, organized as 3+2*n+3 (2 groups (means, LFCs) X n conditions)
                    # ZINB header line has names of conditions, organized as 3+4*n+3 (4 groups X n conditions)
                    if self.inputs.filetype == "anova":
                        n = int((len(w) - 6) / 2)
                    elif self.inputs.filetype == "zinb":
                        n = int((len(headers) - 6) / 4)
                    headers = headers[3 : 3 + n]
                    headers = [x.replace("Mean_", "") for x in headers]
                else:
                    means = [
                        float(x) for x in w[3 : 3 + n]
                    ]  # take just the columns of means
                    lfcs = [
                        float(x) for x in w[3 + n : 3 + n + n]
                    ]  # take just the columns of LFCs
                    qval = float(w[-2])
                    data.append((w, means, lfcs, qval))

        data.sort(key=lambda x: x[-1])
        hits, LFCs = [], []
        for k, (w, means, lfcs, qval) in enumerate(data):
            if (self.inputs.top_k == -1 and qval < self.inputs.adj_p_value) or (
                self.inputs.top_k != -1 and k < self.inputs.top_k
            ):
                mm = round(numpy.mean(means), 1)
                if mm < self.inputs.low_mean_filter:
                    print("excluding %s/%s, mean(means)=%s" % (w[0], w[1], mm))
                else:
                    hits.append(w)
                    LFCs.append(lfcs)

        print("heatmap based on %s genes" % len(hits))
        gene_names = ["%s/%s" % (w[0], w[1]) for w in hits]
        hash = {}
        headers = [h.replace("Mean_", "") for h in headers]
        for i, col in enumerate(headers):
            hash[col] = FloatVector([x[i] for x in LFCs])
        df = DataFrame(hash)
        heatmapFunc = self.make_heatmapFunc()
        heatmapFunc(df, StrVector(gene_names), self.inputs.output_path)

    def make_heatmapFunc(self):
        r(
            """
make_heatmap = function(lfcs,gene_names,outfilename) { 
rownames(lfcs) = gene_names
suppressMessages(require(gplots))
colors <- colorRampPalette(c("red", "white", "blue"))(n = 200)

C = length(colnames(lfcs))
R = length(rownames(lfcs))
W = 300+C*30
H = 300+R*15

png(outfilename,width=W,height=H)
#defaults are lwid=lhei=c(1.5,4)
#heatmap.2(as.matrix(lfcs),col=colors,margin=c(12,12),lwid=c(2,6),lhei=c(0.1,2),trace="none",cexCol=1.4,cexRow=1.4,key=T) # make sure white=0
#heatmap.2(as.matrix(lfcs),col=colors,margin=c(12,12),trace="none",cexCol=1.2,cexRow=1.2,key=T) # make sure white=0 # setting margins was causing failures, so remove it 8/22/21
heatmap.2(as.matrix(lfcs),col=colors,margin=c(12,12),trace="none",cexCol=1.2,cexRow=1.2,key=T) # actually, margins was OK, so the problem must have been with lhei and lwid
dev.off()
}
      """
        )
        return globalenv["make_heatmap"]

    


if __name__ == "__main__":

    (args, kwargs) = transit_tools.clean_args(sys.argv[1:])

    G = Norm.from_args(sys.argv[1:])
    G.Run()
