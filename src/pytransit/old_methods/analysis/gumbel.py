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
import warnings

from pytransit.old_methods import analysis_base as base
from pytransit.tools import transit_tools
from pytransit.tools import tnseq_tools
from pytransit.tools import norm_tools
from pytransit.tools import stat_tools


############# GUI ELEMENTS ##################

short_name = "gumbel"
long_name = "Gumbel"
short_desc = "Bayesian analysis of essentiality based on long gaps."
long_desc = """Bayesian methods of analyzing longest runs of non-insertions in a row. Estimates the parameters using the MCMC sampling, and estimates posterior probabilities of essentiality. 

Reference: DeJesus et al. (2013; Bioinformatics)"""
transposons = ["himar1"]
columns = ["Orf", "Name", "Desc", "k", "n", "r", "s", "zbar", "Call"]


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
            GumbelMethod,
            GumbelGUI,
            [GumbelFile],
        )


################## FILE ###################


class GumbelFile(base.TransitFile):
    def __init__(self):
        base.TransitFile.__init__(self, "#Gumbel", columns)

    def get_header(self, path):
        ess = 0
        unc = 0
        non = 0
        short = 0
        with open(path) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                tmp = line.strip().split("\t")
                if tmp[-1] == "E":
                    ess += 1
                if tmp[-1] == "U":
                    unc += 1
                if tmp[-1] == "NE":
                    non += 1
                if tmp[-1] == "S":
                    short += 1

        text = """Results:
    Essentials: %s
    Uncertain: %s
    Non-Essential: %s
    Short: %s
        """ % (
            ess,
            unc,
            non,
            short,
        )
        return text


################# GUI ##################


class GumbelGUI(base.AnalysisGUI):
    def define_panel(self, wxobj):
        self.wxobj = wxobj
        gumbelPanel = wx.Panel(
            self.wxobj,
            wx.ID_ANY,
            wx.DefaultPosition,
            wx.DefaultSize,
            wx.TAB_TRAVERSAL,
        )

        gumbelSection = wx.BoxSizer(wx.VERTICAL)

        gumbelLabel = wx.StaticText(
            gumbelPanel, wx.ID_ANY, "Gumbel Options", wx.DefaultPosition, (130, 20)
        )

        gumbelLabel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        gumbelSection.Add(gumbelLabel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        mainSizer1 = wx.BoxSizer(wx.VERTICAL)

        # Samples
        (
            gumbelSampleLabel,
            self.wxobj.gumbelSampleText,
            sampleSizer,
        ) = self.defineTextBox(
            gumbelPanel,
            "Samples:",
            "10000",
            "These are the number of samples to take when estimating the parameters. More samples give more accurate estimates of the parameters at the cost of computation time.",
        )
        mainSizer1.Add(sampleSizer, 1, wx.EXPAND, 5)

        # Burn-In
        (
            gumbelBurninLabel,
            self.wxobj.gumbelBurninText,
            burninSizer,
        ) = self.defineTextBox(
            gumbelPanel,
            "Burn-In:",
            "500",
            "These are the number of samples to take before  beginning to estimate the parameters. Allows the MCMC sampler to 'converge' to the true parameter space. More samples give more accurate estimates of the parameters at the cost of computation time.",
        )
        mainSizer1.Add(burninSizer, 1, wx.EXPAND, 5)

        # Trim
        (gumbelTrimLabel, self.wxobj.gumbelTrimText, trimSizer) = self.defineTextBox(
            gumbelPanel,
            "Trim:",
            "1",
            "The MCMC sample will keep every i-th sample. A value of '1' will take all samples. Larger values will reduces autocorrelation at the cost of a substantial cost in computation time.",
        )
        mainSizer1.Add(trimSizer, 1, wx.EXPAND, 5)

        # Min Read
        gumbelReadChoiceChoices = ["1", "2", "3", "4", "5"]
        (
            gumbelReadLabel,
            self.wxobj.gumbelReadChoice,
            readSizer,
        ) = self.defineChoiceBox(
            gumbelPanel,
            "Minimum Read:",
            gumbelReadChoiceChoices,
            "This is the minimum number of reads to consider a 'true' insertion. Value of 1 will consider all insertions. Larger values allow the method to ignore spurious insertions which might interrupt a run of non-insertions. Noisy datasets or those with many replicates can beneffit from increasing this.",
        )
        mainSizer1.Add(readSizer, 1, wx.EXPAND, 5)

        # Replicates
        gumbelRepChoiceChoices = ["Sum", "Mean"]
        (gumbelRepLabel, self.wxobj.gumbelRepChoice, repSizer) = self.defineChoiceBox(
            gumbelPanel,
            "Replicates:",
            gumbelRepChoiceChoices,
            "Determines how to handle replicates, and their read-counts. When using many replicates, summing read-counts may make spurious counts appear to be significantly large and interrupt a run of non-insertions.",
        )
        mainSizer1.Add(repSizer, 1, wx.EXPAND, 5)

        gumbelSection.Add(mainSizer1, 1, wx.EXPAND, 5)

        gumbelButton = wx.Button(
            gumbelPanel, wx.ID_ANY, "Run Gumbel", wx.DefaultPosition, wx.DefaultSize, 0
        )
        gumbelSection.Add(gumbelButton, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        gumbelPanel.SetSizer(gumbelSection)
        gumbelPanel.Layout()
        gumbelSection.Fit(gumbelPanel)

        # Connect events
        gumbelButton.Bind(wx.EVT_BUTTON, self.wxobj.RunMethod)

        self.panel = gumbelPanel


########## METHOD #######################

EXACT = 20

ALPHA = 1
BETA = 1


class GumbelMethod(base.SingleConditionMethod):
    usage_string = """python3 %s gumbel <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -b <integer>    :=  Number of Burn-in samples. Default -b 500
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -t <integer>    :=  Trims all but every t-th value. Default: -t 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        -iN <float>     :=  Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 0
        """ % sys.argv[0]
    
    def __init__(
        self,
        ctrldata,
        annotation_path,
        output_file,
        samples=10000,
        burnin=500,
        trim=1,
        minread=1,
        replicates="Sum",
        normalization=None,
        LOESS=False,
        ignore_codon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
    ):

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
        self.samples = samples
        self.burnin = burnin
        self.trim = trim
        self.minread = minread

        self.cache_nn = {}
        self.cache_expruns = {}

    @classmethod
    def from_gui(self, wxobj):
        """ """

        # Get Annotation file
        from pytransit.universal_data import universal
        annotation_path = universal.annotation_path
        if not transit_tools.validate_annotation(annotation_path):
            return None

        # Get selected files
        ctrldata = wxobj.ctrlSelected()
        if not transit_tools.validate_control_datasets(ctrldata):
            return None

        # Validate transposon types
        if not transit_tools.validate_transposons_used(ctrldata, transposons):
            return None

        # Read the parameters from the wxPython widgets
        try:
            minread = int(
                wxobj.gumbelReadChoice.GetString(
                    wxobj.gumbelReadChoice.GetCurrentSelection()
                )
            )
        except:
            warnings.warn(
                "Warning: problem reading minimum read parameter. Assuming a value of '1'"
            )
            minread = 1

        try:
            samples = int(wxobj.gumbelSampleText.GetValue())
        except:
            warnings.warn(
                "Warning: problem reading samples parameter. Assuming a value of '10000'"
            )
            samples = 10000

        try:
            burnin = int(wxobj.gumbelBurninText.GetValue())
        except:
            warnings.warn(
                "Warning: problem reading burnin parameter. Assuming a value of '500'"
            )
            burnin = 500

        try:
            trim = int(wxobj.gumbelTrimText.GetValue())
        except:
            warnings.warn(
                "Warning: problem reading trim parameter. Assuming a value of '1'"
            )
            trim = 1

        try:
            replicates = wxobj.gumbelRepChoice.GetString(
                wxobj.gumbelRepChoice.GetCurrentSelection()
            )
        except:
            warnings.warn(
                "Warning: problem reading replicates parameter. Assuming a value of 'Mean'"
            )
            replicates = "Mean"

        ignore_codon = True

        try:
            n_terminus = float(wxobj.globalNTerminusText.GetValue())
        except:
            warnings.warn(
                "Warning: problem reading n_terminus parameter. Assuming a value of '0.00'"
            )
            n_terminus = 0.0

        try:
            c_terminus = float(wxobj.globalCTerminusText.GetValue())
        except:
            warnings.warn(
                "Warning: problem reading c_terminus parameter. Assuming a value of '0.00'"
            )
            c_terminus = 0.0

        normalization = None
        LOESS = False

        # Get output path
        name = transit_tools.basename(ctrldata[0])
        defaultFileName = "gumbel_%s_s%d_b%d_t%d.dat" % (
            ".".join(name.split(".")[:-1]),
            samples,
            burnin,
            trim,
        )
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path:
            return None
        output_file = open(output_path, "w")

        return self(
            ctrldata,
            annotation_path,
            output_file,
            samples,
            burnin,
            trim,
            minread,
            replicates,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
            wxobj,
        )

    @classmethod
    def from_args(self, args, kwargs):

        ctrldata = args[0].split(",")
        annotation_path = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        samples = int(kwargs.get("s", 10000))
        burnin = int(kwargs.get("b", 500))
        trim = int(kwargs.get("t", 1))
        minread = int(kwargs.get("m", 1))
        replicates = kwargs.get("r", "Sum")
        normalization = None
        LOESS = False
        ignore_codon = True
        n_terminus = float(kwargs.get("iN", 0.0))
        c_terminus = float(kwargs.get("iC", 0.0))

        return self(
            ctrldata,
            annotation_path,
            output_file,
            samples,
            burnin,
            trim,
            minread,
            replicates,
            normalization,
            LOESS,
            ignore_codon,
            n_terminus,
            c_terminus,
        )

    def Run(self):

        logging.log("Starting Gumbel Method")

        # Set Default parameter values
        w1 = 0.15
        w0 = 1.0 - w1
        ALPHA = 1
        BETA = 1
        ALPHA_w = 600
        BETA_w = 3400
        mu_c = 0
        acctot = 0.0
        phi_start = 0.3
        sigma_c = 0.01

        start_time = time.time()

        

        # Get orf data
        logging.log("Reading Annotation")

        # Validate data has empty sites
        # (status, genome) = transit_tools.validate_wig_format(self.ctrldata, wxobj=self.wxobj)
        # if status <2: tn_used = "himar1"
        # else: tn_used = "tn5"

        logging.log("Getting Data")
        (data, position) = transit_tools.get_validated_data(
            self.ctrldata, wxobj=self.wxobj
        )
        (K, N) = data.shape
        merged = numpy.sum(data, axis=0)
        nsites, nzeros = (
            merged.shape[0],
            numpy.sum(merged == 0),
        )  # perhaps I should say >minCount
        sat = (nsites - nzeros) / float(nsites)

        if self.normalization and self.normalization != "nonorm":
            logging.log("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(
                data, self.normalization, self.ctrldata, self.annotation_path
            )

        G = tnseq_tools.Genes(
            self.ctrldata,
            self.annotation_path,
            minread=self.minread,
            reps=self.replicates,
            ignore_codon=self.ignore_codon,
            n_terminus=self.n_terminus,
            c_terminus=self.c_terminus,
            data=data,
            position=position,
        )

        ii_good = numpy.array(
            [self.good_orf(g) for g in G]
        )  # Gets index of the genes that can be analyzed

        K = G.local_insertions()[ii_good]
        N = G.local_sites()[ii_good]
        R = G.local_runs()[ii_good]
        S = G.local_gap_span()[ii_good]
        T = G.local_gene_span()[ii_good]

        logging.log("Doing Regression")
        mu_s, temp, sigma_s = stat_tools.regress(
            R, S
        )  # Linear regression to estimate mu_s, sigma_s for span data
        mu_r, temp, sigma_r = stat_tools.regress(
            S, R
        )  # Linear regression to estimate mu_r, sigma_r for run data

        N_GENES = len(G)
        N_GOOD = sum(ii_good)

        logging.log("Setting Initial Class")
        Z_sample = numpy.zeros((N_GOOD, self.samples))
        Z = [self.classify(g.n, g.r, 0.5) for g in G if self.good_orf(g)]
        Z_sample[:, 0] = Z
        N_ESS = numpy.sum(Z_sample[:, 0] == 1)

        phi_sample = numpy.zeros(self.samples)  # []
        phi_sample[0] = phi_start
        phi_old = phi_start
        phi_new = 0.00

        SIG = numpy.array(
            [
                self.sigmoid(g.s, g.t) * scipy.stats.norm.pdf(g.r, mu_r * g.s, sigma_r)
                for g in G
                if self.good_orf(g)
            ]
        )

        i = 1
        count = 0
        while i < self.samples:

            try:
                # PHI
                acc = 1.0
                phi_new = phi_old + random.gauss(mu_c, sigma_c)
                i0 = Z_sample[:, i - 1] == 0
                if (
                    phi_new > 1
                    or phi_new <= 0
                    or (
                        self.f_non(phi_new, N[i0], R[i0])
                        - self.f_non(phi_old, N[i0], R[i0])
                    )
                    < math.log(random.uniform(0, 1))
                ):
                    phi_new = phi_old
                    acc = 0.0
                    flag = 0

                # Z
                Z = self.sample_Z(phi_new, w1, N, R, S, T, mu_s, sigma_s, SIG)

                # w1
                N_ESS = sum(Z == 1)
                w1 = scipy.stats.beta.rvs(N_ESS + ALPHA_w, N_GOOD - N_ESS + BETA_w)

                count += 1
                acctot += acc

                if (count > self.burnin) and (count % self.trim == 0):
                    phi_sample[i] = phi_new
                    Z_sample[:, i] = Z
                    i += 1

            except ValueError as e:
                logging.log("Error: %s" % e)
                logging.log("This is likely to have been caused by poor data (e.g. too sparse)." )
                logging.log("If the density of the dataset is too low, the Gumbel method will not work.")
                logging.log("Quitting.")
                return

            #            print(i,phi_new,w1,G[idxG].name,N[idxN],R[idxN],Z[idxN])

            phi_old = phi_new
            # Update progress
            percentage = (100.0 * (count + 1) / (self.samples + self.burnin))
            text = "Running Gumbel Method with Binomial Essentiality Calls... %5.1f%%" % percentage
            progress_update(text, percentage)

        ZBAR = numpy.apply_along_axis(numpy.mean, 1, Z_sample)
        (ess_t, non_t) = stat_tools.bayesian_essentiality_thresholds(ZBAR)
        binomial_n = math.log10(0.05) / math.log10(G.global_phi())

        # Orf    k   n   r   s   zbar
        self.output.write("#Gumbel\n")
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
                "#GUI with: ctrldata=%s, annotation=%s, output=%s, samples=%s, minread=%s, trim=%s\n"
                % (
                    ",".join(self.ctrldata).encode("utf-8"),
                    self.annotation_path.encode("utf-8"),
                    self.output.name.encode("utf-8"),
                    self.samples,
                    self.minread,
                    self.trim,
                )
            )
        else:
            self.output.write("#Console: python3 %s\n" % " ".join(sys.argv))

        self.output.write("#Data: %s\n" % (",".join(self.ctrldata).encode("utf-8")))
        self.output.write("#Annotation path: %s\n" % self.annotation_path.encode("utf-8") )
        self.output.write("#Trimming of TAs near termini: N-term=%s, C-term=%s (fraction of ORF length)\n" % (self.n_terminus, self.c_terminus) )
        self.output.write("#Significance thresholds (FDR-corrected): Essential if Zbar>%f, Non-essential if Zbar<%f\n" % (ess_t, non_t) )
        self.output.write("#Metropolis-Hastings Acceptance-Rate:\t%2.2f%%\n" % (100.0 * acctot / count) )
        self.output.write("#Total Iterations Performed:\t%d\n" % count)
        self.output.write("#Sample Size:\t%d\n" % i)
        self.output.write("#Total number of TA sites: %s\n" % nsites)
        self.output.write("#Genome-wide saturation: %s\n" % (round(sat, 3)) )  # datasets merged
        self.output.write("#phi estimate:\t%f (non-insertion probability in non-essential regions)\n" % numpy.average(phi_sample) )
        self.output.write("#Minimum number of TA sites with 0 insertions to be classified as essential by Binomial: \t%0.3f\n" % binomial_n )
        self.output.write("#Time: %s s\n" % (round(time.time() - start_time, 1)))

        i = 0
        data, calls = [], []
        for j, g in enumerate(G):
            if not self.good_orf(g):
                zbar = -1.0
            else:
                zbar = ZBAR[i]
                i += 1
            if zbar > ess_t:
                call = "E"
            elif G.local_sites()[j] > binomial_n and G.local_thetas()[j] == 0.0:
                call = "EB"
            elif non_t <= zbar <= ess_t:
                call = "U"
            elif 0 <= zbar < non_t:
                call = "NE"
            else:
                call = "S"

            data.append(
                "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%s\n"
                % (g.orf, g.name, g.desc, g.k, g.n, g.r, g.s, zbar, call)
            )
            calls.append(call)
        data.sort()

        self.output.write("#Summary of Essentiality Calls:\n")
        self.output.write("#  E  = %4s (essential based on Gumbel)\n" % (calls.count("E")))
        self.output.write("#  EB = %4s (essential based on Binomial)\n" % (calls.count("EB")))
        self.output.write("#  NE = %4s (non-essential)\n" % (calls.count("NE")))
        self.output.write("#  U  = %4s (uncertain)\n" % (calls.count("U")))
        self.output.write("#  S  = %4s (too short)\n" % (calls.count("S")))

        self.output.write("#%s\n" % "\t".join(columns))
        for line in data:
            self.output.write(line)
        self.output.close()
        logging.log("")  # Printing empty line to flush stdout
        logging.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)
        self.finish()
        logging.log("Finished Gumbel Method")

    def good_orf(self, gene):
        return gene.n >= 3 and gene.t >= 150

    def expected_runs_cached(self, n, q):
        if (n, q) not in self.cache_expruns:
            self.cache_expruns[(n, q)] = tnseq_tools.expected_runs(n, q)
        return self.cache_expruns[(n, q)]

    def classify(self, n, r, p):
        if n == 0:
            return 0
        q = 1 - p
        B = 1 / math.log(1 / p)
        u = math.log(n * q, 1 / p)
        BetaGamma = B * tnseq_tools.get_gamma()
        if (
            n < EXACT
        ):  # estimate more accurately based on expected run len, using exact calc for small genes
            exprun = self.expected_runs_cached(n, p)
            u = (
                exprun - BetaGamma
            )  # u is mu of Gumbel (mean=mu+gamma*beta); matching of moments
            # https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/
        pval = 1 - numpy.exp(scipy.stats.gumbel_r.logcdf(r, u, B))
        if pval < 0.05:
            return 1
        else:
            return 0

    def f_non(self, p, N, R):  # pass in P_nonins as p
        q = 1.0 - p
        BetaGamma = tnseq_tools.get_gamma() / math.log(1 / p)
        total = numpy.log(scipy.stats.beta.pdf(p, ALPHA, BETA))
        mu = numpy.log(N * q) / numpy.log(1 / p)
        for i in range(
            len(N)
        ):  # estimate more accurately based on expected run len, using exact calc for small genes
            if N[i] < EXACT:
                mu[i] = self.expected_runs_cached(int(N[i]), p) - BetaGamma
        sigma = 1 / math.log(1 / p)
        # for i in range(len(N)): print('\t'.join([str(x) for x in N[i],R[i],self.expected_runs_cached(int(N[i]),q),mu[i],scipy.stats.gumbel_r.pdf(R[i], mu[i], sigma)]))
        total += numpy.sum(scipy.stats.gumbel_r.logpdf(R, mu, sigma))
        return total

    def sample_Z(self, p, w1, N, R, S, T, mu_s, sigma_s, SIG):
        G = len(N)
        q = 1.0 - p
        BetaGamma = tnseq_tools.get_gamma() / math.log(1 / p)
        mu = numpy.log(N * q) / numpy.log(1 / p)
        for i in range(
            len(N)
        ):  # estimate more accurately based on expected run len, using exact calc for small genes
            if N[i] < EXACT:
                mu[i] = self.expected_runs_cached(int(N[i]), p) - BetaGamma
        sigma = 1.0 / math.log(1.0 / p)
        h0 = (
            (numpy.exp(scipy.stats.gumbel_r.logpdf(R, mu, sigma)))
            * scipy.stats.norm.pdf(S, mu_s * R, sigma_s)
            * (1 - w1)
        )
        h1 = SIG * w1
        h1 += 1e-10
        h0 += 1e-10  # to prevent div-by-zero; if neither class is probable, p(z1) should be ~0.5
        p_z1 = h1 / (h0 + h1)
        return scipy.stats.binom.rvs(1, p_z1, size=G)

    def sigmoid(self, d, n):
        Kn = 0.1
        MEAN_DOMAIN_SPAN = 300

        if d == 0:
            return 0.00
        f = 1.0 / (1.0 + math.exp(Kn * (MEAN_DOMAIN_SPAN - d)))
        # if n in self.cache_nn: return f/self.cache_nn[n]
        tot = 0
        N = int(n + 1)
        for i in range(1, N):
            tot += 1.0 / (1.0 + math.exp(Kn * (MEAN_DOMAIN_SPAN - i)))
        self.cache_nn[n] = tot
        return f / tot
