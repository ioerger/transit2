import sys
import os
import time
import math
import random
import datetime

import numpy

from pytransit.old_methods import analysis_base as base
from pytransit.tools.transit_tools import HAS_R, r, DataFrame, globalenv, IntVector, FloatVector, StrVector, rpackages
from pytransit.tools import logging, tnseq_tools

############# Description ##################

short_name  = "corrplot"
long_name   = "Corrplot"
short_desc  = "Correlation among TnSeq datasets"
long_desc   = "Correlation among TnSeq datasets"
transposons = ["himar1"  , "tn5"]
columns     = ["Position", "Reads", "Genes"]


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
            CorrplotMethod,
            CorrplotGUI,
            [],
        )  # [CorrplotFile])


################## FILE ###################

# there is no output file that could be loaded into the GUI

# class CorrplotFile(base.TransitFile):
#
#    def __init__(self):
#        base.TransitFile.__init__(self, "#CombinedWig", columns)
#
#    def get_header(self, path):
#        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
#        return text

################# GUI ##################

# right now, tnseq_stats is just intended for the command-line; TRI


class CorrplotGUI(base.AnalysisGUI):
    def __init__(self):
        base.AnalysisGUI.__init__(self)


########## METHOD #######################


class CorrplotMethod(base.SingleConditionMethod):
    """   
    Norm
 
    """

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
        if not HAS_R:
            raise Exception(f'''
                Error: R and rpy2 (~= 3.0) required to run corrplot.
                After installing R, you can install rpy2 using the command \"pip install 'rpy2~=3.0'\"
            ''')
        
        if kwargs.get("-help", False):
            print(self.usage_string)
            sys.exit(0)
        if len(args) < 2:
            print(self.usage_string)
            sys.exit(0)
        self.gene_means = args[0]
        self.outfile = args[1]
        self.filetype = "gene_means"
        if "-anova" in kwargs or "anova" in kwargs:
            self.filetype = "anova"
        if "-zinb" in kwargs or "zinb" in kwargs:
            self.filetype = "zinb"
        return self(self.gene_means, outfile=self.outfile)

    def Run(self):

        logging.log("Starting Corrplot")
        start_time = time.time()

        # assume first non-comment line is header; samples are
        headers = None
        data, means = [], []

        if self.filetype == "gene_means":
            with open(self.gene_means) as file:
                for line in file:
                    w = line.rstrip().split("\t")
                    if line[0] == "#":
                        headers = w[3:]
                        continue  # last comment line has names of samples
                    data.append(w)
                    cnts = [float(x) for x in w[3:]]
                    means.append(cnts)
        elif self.filetype == "anova" or self.filetype == "zinb":
            n = -1  # number of conditions
            with open(self.gene_means) as file:
                for line in file:
                    w = line.rstrip().split("\t")
                    if line[0] == "#" or (
                        "pval" in line and "padj" in line
                    ):  # check for 'pval' for backwards compatibility
                        headers = w
                        continue  # keep last comment line as headers
                    if n == -1:
                        # ANOVA header line has names of conditions, organized as 3+2*n+3 (2 groups (means, LFCs) X n conditions)
                        # ZINB header line has names of conditions, organized as 3+4*n+3 (4 groups X n conditions)
                        if self.filetype == "anova":
                            n = int((len(w) - 6) / 2)
                        elif self.filetype == "zinb":
                            n = int((len(headers) - 6) / 4)
                        headers = headers[3 : 3 + n]
                        headers = [x.replace("Mean_", "") for x in headers]
                    vals = [
                        float(x) for x in w[3 : 3 + n]
                    ]  # take just the columns of means
                    qval = float(w[-2])
                    if qval < 0.05:
                        data.append(w)
                        means.append(vals)
        else:
            print("filetype not recognized: %s" % self.filetype)
            sys.exit(-1)
        print("correlations based on %s genes" % len(means))

        genenames = ["%s/%s" % (w[0], w[1]) for w in data]
        hash = {}
        headers = [h.replace("Mean_", "") for h in headers]
        headers = [h.replace("-", ".") for h in headers]  # because of R conversion
        headers = [
            "X" + x if x[0].isdigit() else x for x in headers
        ]  # because R prepends 'X' to column names starting with a digit
        for i, col in enumerate(headers):
            hash[col] = FloatVector([x[i] for x in means])
        df = DataFrame(hash)  # can't figure out how to set rownames

        corrplotFunc = self.make_corrplotFunc()
        corrplotFunc(
            df, StrVector(headers), StrVector(genenames), self.outfile
        )  # pass headers to put cols in order, since df comes from dict

        self.finish()
        logging.log("Finished Corrplot")

    def make_corrplotFunc(self):
        r(
            """ # R function...
make_corrplot = function(means,headers,genenames,outfilename) { 
means = means[,headers] # put cols in correct order
rownames(means) = genenames
suppressMessages(require(corrplot))
png(outfilename)
corrplot(cor(means))
dev.off()
}
      """
        )
        return globalenv["make_corrplot"]

    usage_string = "usage: python3 %s corrplot <gene_means> <output.png> [-anova|-zinb]" % sys.argv[0]
        # could add a flag for padj cutoff (or top n most signif genes)


if __name__ == "__main__":

    (args, kwargs) = transit_tools.clean_args(sys.argv[1:])

    G = Norm.from_args(sys.argv[1:])
    G.Run()
