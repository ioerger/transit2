from pytransit.components.parameter_panel import panel, progress_update
import pytransit.components.results_area as results_area
import sys

from pytransit.tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub

import os
import time
import math
import random
import numpy
import scipy.stats
import datetime

from pytransit.methods import analysis_base as base
import pytransit.tools.transit_tools as transit_tools
import pytransit.tools.tnseq_tools as tnseq_tools
import pytransit.tools.norm_tools as norm_tools
import pytransit.tools.stat_tools as stat_tools

# method_name = "hmm"


############# GUI ELEMENTS ##################

Analysis.short_name = "hmm"
Analysis.long_name = "HMM"
Analysis.short_desc = "Analysis of genomic regions using a Hidden Markov Model"
Analysis.long_desc = """Analysis of essentiality in the entire genome using a Hidden Markov Model. Capable of determining regions with different levels of essentiality representing Essential, Growth-Defect, Non-Essential and Growth-Advantage regions.

Reference: DeJesus et al. (2013; BMC Bioinformatics)
"""
Analysis.transposons = ["himar1"]
Analysis.columns_sites = [
    "Location",
    "Read Count",
    "Probability - ES",
    "Probability - GD",
    "Probability - NE",
    "Probability - GA",
    "State",
    "Gene",
]
Analysis.columns_genes = [
    "Orf",
    "Name",
    "Description",
    "Total Sites",
    "Num. ES",
    "Num. GD",
    "Num. NE",
    "Num. GA",
    "Avg. Insertions",
    "Avg. Reads",
    "State Call",
]


############# Analysis Method ##############


class Analysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(
            self,
            Analysis.short_name,
            Analysis.long_name,
            Analysis.short_desc,
            Analysis.long_desc,
            Analysis.transposons,
            HMMMethod,
            HMMGUI,
            [HMMSitesFile, HMMGenesFile],
        )


################## FILE ###################


class HMMSitesFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#HMM - Sites", Analysis.columns_sites)

    def get_header(self, path):
        es = 0
        gd = 0
        ne = 0
        ga = 0
        T = 0
        with open(path) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                tmp = line.strip().split("\t")
                if len(tmp) == 7:
                    col = -1
                else:
                    col = -2
                if tmp[col] == "ES":
                    es += 1
                elif tmp[col] == "GD":
                    gd += 1
                elif tmp[col] == "NE":
                    ne += 1
                elif tmp[col] == "GA":
                    ga += 1
                else:
                    print(tmp)
                T += 1

        text = """Results:
    Essential: %1.1f%%
    Growth-Defect: %1.1f%%
    Non-Essential: %1.1f%%
    Growth-Advantage: %1.1f%%
            """ % (
            100.0 * es / T,
            100.0 * gd / T,
            100.0 * ne / T,
            100.0 * ga / T,
        )
        return text


class HMMGenesFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#HMM - Genes", Analysis.columns_genes)

    def get_header(self, path):
        es = 0
        gd = 0
        ne = 0
        ga = 0
        T = 0
        with open(path) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                tmp = line.strip().split("\t")
                if len(tmp) < 5:
                    continue
                if tmp[-1] == "ES":
                    es += 1
                if tmp[-1] == "GD":
                    gd += 1
                if tmp[-1] == "NE":
                    ne += 1
                if tmp[-1] == "GA":
                    ga += 1

        text = """Results:
    Essential: %s
    Growth-Defect: %s
    Non-Essential: %s
    Growth-Advantage: %s
            """ % (
            es,
            gd,
            ne,
            ga,
        )

        return text


############# GUI ##################


class HMMGUI(base.AnalysisGUI):
    pass


########## CLASS #######################


class HMMMethod(base.SingleConditionMethod):
    

    def __init__(
        self,
        ctrldata,
        annotation_path,
        output_file,
        replicates="Mean",
        normalization=None,
        loess=False,
        ignore_codon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
    ):

        base.SingleConditionMethod.__init__(
            self,
            Analysis.short_name,
            Analysis.long_name,
            Analysis.short_desc,
            Analysis.long_desc,
            ctrldata,
            annotation_path,
            output_file,
            replicates=replicates,
            normalization=normalization,
            loess=loess,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            wxobj=wxobj,
        )

        try:
            with open(ctrldata[0]) as file:
                T = len(
                    [
                        1
                        for line in file.readlines()
                        if not line.startswith("#")
                    ]
                )
                self.maxiterations = T * 4 + 1
        except:
            self.maxiterations = 100
        self.count = 1

    @classmethod
    def from_gui(self, wxobj):
        pass

    @classmethod
    def from_args(self, rawargs):
        pass

    def Run(self):
        # Get data
        transit_tools.log("Getting Data")
        (data, position) = transit_tools.get_validated_data(
            self.inputs.ctrldata, wxobj=self.wxobj
        )
        (K, N) = data.shape

        # Normalize data
        if self.inputs.normalization != "nonorm":
            transit_tools.log("Normalizing using: %s" % self.inputs.normalization)
            wig_list = [] # FIXME: normalize_data() needs to be rewritten to not need a wig_list
            (data, factors) = norm_tools.normalize_data(
                data, self.inputs.normalization, annotationPath=self.inputs.annotation_path
            )

        # Do loess
        if self.inputs.loess:
            transit_tools.log("Performing loess Correction")
            for j in range(K):
                data[j] = stat_tools.loess_correction(position, data[j])

        hash = transit_tools.get_pos_hash(self.inputs.annotation_path)
        rv2info = transit_tools.get_gene_info(self.inputs.annotation_path)

        if len(self.inputs.ctrldata) > 1:
            transit_tools.log("Combining Replicates as '%s'" % self.inputs.replicates)
        O = (
            tnseq_tools.combine_replicates(data, method=self.inputs.replicates) + 1
        )  # Adding 1 to because of shifted geometric in scipy

        # Parameters
        Nstates = 4
        label = {0: "ES", 1: "GD", 2: "NE", 3: "GA"}

        reads = O - 1
        reads_nz = sorted(reads[reads != 0])
        size = len(reads_nz)
        mean_r = numpy.average(reads_nz[: int(0.95 * size)])
        mu = numpy.array([1 / 0.99, 0.01 * mean_r + 2, mean_r, mean_r * 5.0])
        # mu = numpy.array([1/0.99, 0.1 * mean_r + 2,  mean_r, mean_r*5.0])
        L = 1.0 / mu
        B = []  # Emission Probability Distributions
        for i in range(Nstates):
            B.append(scipy.stats.geom(L[i]).pmf)

        pins = self.calculate_pins(O - 1)
        pins_obs = sum([1 for rd in O if rd >= 2]) / float(len(O))
        pnon = 1.0 - pins
        pnon_obs = 1.0 - pins_obs

        for r in range(100):
            if pnon ** r < 0.01:
                break

        A = numpy.zeros((Nstates, Nstates))
        a = math.log1p(-B[int(Nstates / 2)](1) ** r)
        b = r * math.log(B[int(Nstates / 2)](1)) + math.log(
            1.0 / 3
        )  # change to Nstates-1?
        for i in range(Nstates):
            A[i] = [b] * Nstates
            A[i][i] = a

        PI = numpy.zeros(Nstates)  # Initial state distribution
        PI[0] = 0.7
        PI[1:] = 0.3 / (Nstates - 1)

        

        ###############
        ### VITERBI ###
        (Q_opt, delta, Q) = self.viterbi(A, B, PI, O)
        ###############

        ##################
        ### ALPHA PASS ###
        (log_Prob_Obs, alpha, C) = self.forward_procedure(numpy.exp(A), B, PI, O)
        ##################

        #################
        ### BETA PASS ###
        beta = self.backward_procedure(numpy.exp(A), B, PI, O, C)
        #################

        T = len(O)
        total = 0
        state2count = dict.fromkeys(range(Nstates), 0)
        for t in range(T):
            state = Q_opt[t]
            state2count[state] += 1
            total += 1

        self.output.write("#HMM - Sites\n")
        self.output.write("# Tn-HMM\n")

        if self.wxobj:
            members = sorted(
                [
                    attr
                    for attr in dir(self)
                    if not callable(getattr(self, attr)) and not attr.startswith("__")
                ]
            )
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write(
                "#GUI with: ctrldata=%s, annotation=%s, output=%s\n"
                % (
                    ",".join(self.inputs.ctrldata).encode("utf-8"),
                    self.inputs.annotation_path.encode("utf-8"),
                    self.output.name.encode("utf-8"),
                )
            )
        else:
            self.output.write("#Console: python3 %s\n" % " ".join(sys.argv))

        self.output.write("# \n")
        self.output.write("# Mean:\t%2.2f\n" % (numpy.average(reads_nz)))
        self.output.write("# Median:\t%2.2f\n" % numpy.median(reads_nz))
        self.output.write("# Normalization:\t%s\n" % self.inputs.normalization)
        self.output.write("# loess Correction:\t%s\n" % str(self.inputs.loess))
        self.output.write("# pins (obs):\t%f\n" % pins_obs)
        self.output.write("# pins (est):\t%f\n" % pins)
        self.output.write("# Run length (r):\t%d\n" % r)
        self.output.write("# State means:\n")
        self.output.write(
            "#    %s\n"
            % "   ".join(["%s: %8.4f" % (label[i], mu[i]) for i in range(Nstates)])
        )
        self.output.write("# Self-Transition Prob:\n")
        self.output.write(
            "#    %s\n"
            % "   ".join(["%s: %2.4e" % (label[i], A[i][i]) for i in range(Nstates)])
        )
        self.output.write("# State Emission Parameters (theta):\n")
        self.output.write(
            "#    %s\n"
            % "   ".join(["%s: %1.4f" % (label[i], L[i]) for i in range(Nstates)])
        )
        self.output.write("# State Distributions:")
        self.output.write(
            "#    %s\n"
            % "   ".join(
                [
                    "%s: %2.2f%%" % (label[i], state2count[i] * 100.0 / total)
                    for i in range(Nstates)
                ]
            )
        )

        states = [int(Q_opt[t]) for t in range(T)]
        last_orf = ""
        for t in range(T):
            s_lab = label.get(states[t], "Unknown State")
            gamma_t = (alpha[:, t] * beta[:, t]) / numpy.sum(alpha[:, t] * beta[:, t])
            genes_at_site = hash.get(position[t], [""])
            genestr = ""
            if not (len(genes_at_site) == 1 and not genes_at_site[0]):
                genestr = ",".join(
                    ["%s_(%s)" % (g, rv2info.get(g, "-")[0]) for g in genes_at_site]
                )

            self.output.write(
                "%s\t%s\t%s\t%s\t%s\n"
                % (
                    int(position[t]),
                    int(O[t]) - 1,
                    "\t".join(["%-9.2e" % g for g in gamma_t]),
                    s_lab,
                    genestr,
                )
            )

        self.output.close()

        transit_tools.log("")  # Printing empty line to flush stdout
        transit_tools.log("Finished HMM - Sites Method")
        transit_tools.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)

        # Gene Files
        transit_tools.log("Creating HMM Genes Level Output")
        genes_path = (
            ".".join(self.output.name.split(".")[:-1])
            + "_genes."
            + self.output.name.split(".")[-1]
        )

        tempObs = numpy.zeros((1, len(O)))
        tempObs[0, :] = O - 1
        self.post_process_genes(tempObs, position, states, genes_path)

        transit_tools.log("Adding File: %s" % (genes_path))
        results_area.add(self.output.name)
        self.finish()
        transit_tools.log("Finished HMM Method")

    

if __name__ == "__main__":

    (args, kwargs) = transit_tools.clean_args(sys.argv)

    G = HMMMethod.from_args(sys.argv[1:])

    G.console_message("Printing the member variables:")
    G.print_members()

    print("")
    print("Running:")

    G.Run()
