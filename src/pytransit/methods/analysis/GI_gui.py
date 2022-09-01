from pytransit.components.parameter_panel import panel, progress_update
import sys
import os
import time
import ntpath
import math
import random
import datetime
import heapq
import collections

import numpy
import scipy
import scipy.stats
import heapq
import math
import statsmodels.stats.multitest
from pytransit.basics.lazy_dict import LazyDict

from pytransit.methods import analysis_base as base
from pytransit.tools.transit_tools import wx, pub, basename, HAS_R, FloatVector, DataFrame, StrVector, EOL
from pytransit.methods import analysis_base as base
import pytransit
import pytransit.tools.gui_tools as gui_tools
import pytransit.components.file_display as file_display
import pytransit.tools.transit_tools as transit_tools
import pytransit.tools.console_tools as console_tools
import pytransit.tools.tnseq_tools as tnseq_tools
import pytransit.tools.norm_tools as norm_tools
import pytransit.tools.stat_tools as stat_tools
import pytransit.basics.csv as csv
import pytransit.components.results_area as results_area
from pytransit.universal_data import universal
from pytransit.components.parameter_panel import panel as parameter_panel
from pytransit.components.parameter_panel import panel, progress_update
from pytransit.components.spreadsheet import SpreadSheet
#nfrom pytransit.components.panel_helpers import make_panel, create_run_button, create_button, create_normalization_input, define_choice_box
from pytransit.components.panel_helpers import *

command_name = sys.argv[0]

class Analysis:
    identifier  = "#GI_gui"
    short_name  = "GI_gui"
    long_name   = "GI_gui"
    short_desc  = "Genetic Interaction analysis"
    long_desc   = """Genetic Interaction analysis"""
    transposons = [ "himar1" ]
    
    inputs = LazyDict(
        conditionA1=None,
        conditionB1=None,
        conditionA2=None,
        conditionB2=None,
        combined_wig=None,
        normalization=None,
        output_path=None,
    )
    
    valid_cli_flags = [
        "-n", # normalization flag
        "-o", # output filename (optional)
    ]
    # see gi.py for original usage string and list of flags...
    usage_string = f"""usage: python3 %s GI <combined_wig> <conditionA1> <conditionB1> <conditionA2> <conditionB2> <prot_table> <output_file> [-n <norm>]"""
    
    wxobj = None
    panel = None
    
    def __init__(self, *args, **kwargs):
        Analysis.instance = self
        self.full_name        = f"[{self.short_name}]  -  {self.short_desc}"
        self.transposons_text = transit_tools.get_transposons_text(self.transposons)
        self.filetypes        = [File]
        self.method           = Analysis # backwards compat
        self.gui              = self     # backwards compat
    
    def __str__(self):
        return f"""
            Analysis Method:
                Short Name:  {self.short_name}
                Long Name:   {self.long_name}
                Short Desc:  {self.short_desc}
                Long Desc:   {self.long_desc}
                GUI:         {self.gui}
        """.replace('\n            ','\n').strip()
    
    def __repr__(self):
        return f"{self.inputs}"

    def create_condition_choice(self, panel, sizer, name):
        (
            label,
            ref_condition_wxobj,
            ref_condition_choice_sizer,
        ) = define_choice_box(
            panel,
            name,
            [ "[None]" ] + [x.name for x in universal.session_data.conditions],
            "choose condition",
        )
        sizer.Add(ref_condition_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        return lambda *args: ref_condition_wxobj.GetString(ref_condition_wxobj.GetCurrentSelection())
    
    def create_input_field(self, panel, sizer, label, value,tooltip=None):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text=label,
            default_value=value,
            tooltip_text=tooltip,
        )
        return lambda *args: float(get_text())

    def define_panel(self, _):
        self.panel = make_panel()

        # only need Norm selection and Run button        
        self.value_getters = LazyDict()
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.value_getters.condA1 = self.create_condition_choice(self.panel,main_sizer,"Condition A1:")
        self.value_getters.condB1 = self.create_condition_choice(self.panel,main_sizer,"Condition B1:")
        self.value_getters.condA2 = self.create_condition_choice(self.panel,main_sizer,"Condition A2:")
        self.value_getters.condB2 = self.create_condition_choice(self.panel,main_sizer,"Condition B2:")
        self.value_getters.normalization = create_normalization_input(self.panel, main_sizer) # TTR is default
        self.value_getters.n_terminus = create_n_terminus_input(self.panel, main_sizer)
        self.value_getters.c_terminus = create_c_terminus_input(self.panel, main_sizer)
        self.value_getters.samples = self.create_input_field(self.panel, main_sizer,"Number of samples",10000,"random trials in Monte Carlo simulation")
        self.value_getters.rope = self.create_input_field(self.panel, main_sizer,"ROPE",0.5,"Region of probable equivalence around 0")
        # to do:
        # checkbox [off] for 'correct for genome position bias' (view LOESS fit)
        # checkbox [on] for 'include sites with all zeros'
        # dropdown based on new CL flag for 'significance method': -signif <HDI,prob,BFDR,FWER>
        create_run_button(self.panel, main_sizer)

        parameter_panel.set_panel(self.panel)
        self.panel.SetSizer(main_sizer)
        self.panel.Layout()
        main_sizer.Fit(self.panel)

    @classmethod
    def from_gui(cls, frame):
        with gui_tools.nice_error_log:

            # 
            # get wig files
            # 
            wig_group = universal.session_data.combined_wigs[0] # assume there is only 1 (should check that it has beed defined)
            Analysis.inputs.combined_wig = wig_group.main_path # see components/sample_area.py
            
            # 
            # get annotation
            # 

            Analysis.inputs.annotation_path = universal.session_data.annotation_path
            #if not transit_tools.validate_annotation(Analysis.inputs.annotation): return None 

            
            # 
            # setup custom inputs
            # 
            for each_key, each_getter in Analysis.instance.value_getters.items():
                try:
                    Analysis.inputs[each_key] = each_getter()
                except Exception as error:
                    raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')

            # 
            # save result files
            # 
            Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
                default_file_name="GI_test.dat",
                output_extensions=u'Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
            )
            if not Analysis.inputs.output_path:
                return None

            return Analysis.instance

    @classmethod
    def from_args(cls, args, kwargs):
        console_tools.handle_help_flag(kwargs, cls.usage_string)
        console_tools.handle_unrecognized_flags(cls.valid_cli_flags, kwargs, cls.usage_string)

        normalization = kwargs.get("n", "nonorm") 
        output_path = kwargs.get("o", None)

        # save all the data
        Analysis.inputs.update(dict(
            combined_wig=combined_wig, ###? what if user gives a list of wig files instead of a combined_wig?
            annotation=annotation_path,
            normalization=normalization,
            output_path=output_path,
            n_terminus=0,
            c_terminus=0,
        ))
        
        return Analysis.instance
        
    def Run(self):
        with gui_tools.nice_error_log:
            transit_tools.log("Starting Genetic Interaction analysis")
            start_time = time.time()

            self.inputs.n_terminus = 0
            self.inputs.c_terminus = 0
            self.inputs.samples = 100
            self.inputs.includeZeros = False
            self.inputs.rope=0.5
            self.inputs.signif="HDI"

            # 
            # get data
            # 

            transit_tools.log("Getting Data")
            sites, data, filenames_in_comb_wig = tnseq_tools.read_combined_wig(self.inputs.combined_wig)
            transit_tools.log(f"Normalizing using: {self.inputs.normalization}")
            data, factors = norm_tools.normalize_data(data, self.inputs.normalization)
                
            # # Do LOESS correction if specified
            # if self.LOESS:
            #   transit_tools.log("Performing LOESS Correction")
            #   for j in range(K):
            #     data[j] = stat_tools.loess_correction(position, data[j])


            # 
            # process data
            # 
            # note: self.inputs.condA1 = self.value_getters.condA1()

            transit_tools.log("processing data")
            # get 4 lists of indexes into data (extract columns for 4 conds in comwig)
            # could also pull this info from universal.session_data.samples?
            from pytransit.components.samples_area import sample_table
            indexes = {}
            for i,row in enumerate(sample_table.rows): # defined in components/generic/table.py
              cond = row["condition"]
              if cond not in indexes: indexes[cond] = []
              indexes[cond].append(i)
            condA1 = self.inputs.condA1
            condA2 = self.inputs.condA2
            condB1 = self.inputs.condB1
            condB2 = self.inputs.condB2
            if condA1 not in indexes or len(indexes[condA1])==0: print("error: no samples found for condition %s" % condA1); sys.exit(0)
            if condA2 not in indexes or len(indexes[condA2])==0: print("error: no samples found for condition %s" % condA2); sys.exit(0)
            if condB1 not in indexes or len(indexes[condB1])==0: print("error: no samples found for condition %s" % condB1); sys.exit(0)
            if condB2 not in indexes or len(indexes[condB2])==0: print("error: no samples found for condition %s" % condB2); sys.exit(0)
            print("condA1=%s, indexes=%s" % (condA1,','.join([str(x) for x in indexes[condA1]])))
            print("condA2=%s, indexes=%s" % (condA2,','.join([str(x) for x in indexes[condA2]])))
            print("condB1=%s, indexes=%s" % (condB1,','.join([str(x) for x in indexes[condB1]])))
            print("condB2=%s, indexes=%s" % (condB2,','.join([str(x) for x in indexes[condB2]])))

            dataA1 = data[indexes[condA1]] # select datasets (rows)
            dataA2 = data[indexes[condA2]]
            dataB1 = data[indexes[condB1]]
            dataB2 = data[indexes[condB2]]

            # results: 1 row for each gene; adjusted_label - just a string that gets printed in the header
            (results,adjusted_label) = self.calc_GI(dataA1,dataA2,dataB1,dataB2,sites)

            # 
            # write output
            # 
            # note: first comment line is filetype, last comment line is column headers

            transit_tools.log(f"Adding File: {self.inputs.output_path}")
            results_area.add(self.inputs.output_path)

            # will open and close file, or print to console
            self.print_GI(results,adjusted_label,condA1,condA2,condB1,condB2,sample_table) 

            transit_tools.log("Finished Genetic Interaction analysis")
            transit_tools.log("Time: %0.1fs\n" % (time.time() - start_time))

    # is self.annotation_path already set?

    def calc_GI(self,dataA1,dataA2,dataB1,dataB2,position): # position is vector of TAsite coords

        # Get Gene objects for each condition
        G_A1 = tnseq_tools.Genes(
            [],
            self.inputs.annotation_path,
            data=dataA1, # data[:Na1],
            position=position,
            n_terminus=self.inputs.n_terminus,
            c_terminus=self.inputs.c_terminus,
        )
        G_B1 = tnseq_tools.Genes(
            [],
            self.inputs.annotation_path,
            data=dataB1, # data[Na1 : (Na1 + Nb1)],
            position=position,
            n_terminus=self.inputs.n_terminus,
            c_terminus=self.inputs.c_terminus,
        )
        G_A2 = tnseq_tools.Genes(
            [],
            self.inputs.annotation_path,
            data=dataA2, # data[(Na1 + Nb1) : (Na1 + Nb1 + Na2)],
            position=position,
            n_terminus=self.inputs.n_terminus,
            c_terminus=self.inputs.c_terminus,
        )
        G_B2 = tnseq_tools.Genes(
            [],
            self.inputs.annotation_path,
            data=dataB2, # data[(Na1 + Nb1 + Na2) :],
            position=position,
            n_terminus=self.inputs.n_terminus,
            c_terminus=self.inputs.c_terminus,
        )

        means_list_a1 = []
        means_list_b1 = []
        means_list_a2 = []
        means_list_b2 = []

        var_list_a1 = []
        var_list_a2 = []
        var_list_b1 = []
        var_list_b2 = []

        # Base priors on empirical observations across genes.
        for gene in sorted(G_A1):
            if gene.n > 1:
                A1_data = G_A1[gene.orf].reads.flatten()
                B1_data = G_B1[gene.orf].reads.flatten()
                A2_data = G_A2[gene.orf].reads.flatten()
                B2_data = G_B2[gene.orf].reads.flatten()

                means_list_a1.append(numpy.mean(A1_data))
                var_list_a1.append(numpy.var(A1_data))

                means_list_b1.append(numpy.mean(B1_data))
                var_list_b1.append(numpy.var(B1_data))

                means_list_a2.append(numpy.mean(A2_data))
                var_list_a2.append(numpy.var(A2_data))

                means_list_b2.append(numpy.mean(B2_data))
                var_list_b2.append(numpy.var(B2_data))

        # Priors
        mu0_A1 = scipy.stats.trim_mean(means_list_a1, 0.01)
        mu0_B1 = scipy.stats.trim_mean(means_list_b1, 0.01)
        mu0_A2 = scipy.stats.trim_mean(means_list_a2, 0.01)
        mu0_B2 = scipy.stats.trim_mean(means_list_b2, 0.01)

        s20_A1 = scipy.stats.trim_mean(var_list_a1, 0.01)
        s20_B1 = scipy.stats.trim_mean(var_list_b1, 0.01)
        s20_A2 = scipy.stats.trim_mean(var_list_a2, 0.01)
        s20_B2 = scipy.stats.trim_mean(var_list_b2, 0.01)

        k0 = 1.0
        nu0 = 1.0
        results = [] # results vectors for each gene

        postprob = []
        count = 0
        N = len(G_A1)
        
        # Perform actual analysis
        for gene in G_A1:

            # If there is some data
            if gene.n > 0:
                A1_data = G_A1[gene.orf].reads.flatten()
                B1_data = G_B1[gene.orf].reads.flatten()
                A2_data = G_A2[gene.orf].reads.flatten()
                B2_data = G_B2[gene.orf].reads.flatten()

                #            Time-1   Time-2
                #
                #  Strain-A     A       C
                #
                #  Strain-B     B       D

                try:
                    muA1_post, varA1_post = stat_tools.sample_trunc_norm_post(
                        A1_data, self.inputs.samples, mu0_A1, s20_A1, k0, nu0
                    )
                    muB1_post, varB1_post = stat_tools.sample_trunc_norm_post(
                        B1_data, self.inputs.samples, mu0_B1, s20_B1, k0, nu0
                    )
                    muA2_post, varA2_post = stat_tools.sample_trunc_norm_post(
                        A2_data, self.inputs.samples, mu0_A2, s20_A2, k0, nu0
                    )
                    muB2_post, varB2_post = stat_tools.sample_trunc_norm_post(
                        B2_data, self.inputs.samples, mu0_B2, s20_B2, k0, nu0
                    )

                except Exception as e:
                    muA1_post = varA1_post = numpy.ones(self.inputs.samples)
                    muB1_post = varB1_post = numpy.ones(self.inputs.samples)
                    muA2_post = varA2_post = numpy.ones(self.inputs.samples)
                    muB2_post = varB2_post = numpy.ones(self.inputs.samples)

                logFC_A_post = numpy.log2(muA2_post / muA1_post)
                logFC_B_post = numpy.log2(muB2_post / muB1_post)
                delta_logFC_post = logFC_B_post - logFC_A_post

                alpha = 0.05

                # Get Bounds of the HDI
                l_logFC_A, u_logFC_A = stat_tools.hdi_from_mcmc(logFC_A_post, 1 - alpha)

                l_logFC_B, u_logFC_B = stat_tools.hdi_from_mcmc(logFC_B_post, 1 - alpha)

                l_delta_logFC, u_delta_logFC = stat_tools.hdi_from_mcmc(
                    delta_logFC_post, 1 - alpha
                )

                mean_logFC_A = numpy.mean(logFC_A_post)
                mean_logFC_B = numpy.mean(logFC_B_post)
                mean_delta_logFC = numpy.mean(delta_logFC_post)

                # Is HDI significantly different than ROPE? (i.e. no overlap)
                not_HDI_overlap_bit = (
                    l_delta_logFC > self.inputs.rope or u_delta_logFC < -self.inputs.rope
                )

                # Probability of posterior overlaping with ROPE
                probROPE = numpy.mean(
                    numpy.logical_and(
                        delta_logFC_post >= 0.0 - self.inputs.rope,
                        delta_logFC_post <= 0.0 + self.inputs.rope,
                    )
                )

            # If there is no data, assume empty defaults
            else:
                A1_data = [0, 0]
                B1_data = [0, 0]
                A2_data = [0, 0]
                B2_data = [0, 0]
                muA1_post = varA1_post = numpy.ones(self.inputs.samples)
                muB1_post = varB1_post = numpy.ones(self.inputs.samples)
                muA2_post = varA2_post = numpy.ones(self.inputs.samples)
                muB2_post = varB2_post = numpy.ones(self.inputs.samples)
                logFC_A_post = numpy.log2(muA2_post / muA1_post)
                logFC_B_post = numpy.log2(muB2_post / muB1_post)
                delta_logFC_post = logFC_B_post - logFC_A_post

                mean_logFC_A = 0
                mean_logFC_B = 0
                mean_delta_logFC = 0
                l_logFC_A = 0
                u_logFC_A = 0
                l_logFC_B = 0
                u_logFC_B = 0
                l_delta_logFC = 0
                u_delta_logFC = 0
                probROPE = 1.0

            if numpy.isnan(l_logFC_A):
                l_logFC_A = -10
                u_logFC_A = 10
            if numpy.isnan(l_logFC_B):
                l_logFC_B = -10
                u_logFC_B = 10
            if numpy.isnan(l_delta_logFC):
                l_delta_logFC = -10
                u_delta_logFC = 10

            postprob.append(probROPE)
            results.append(
                (
                    gene.orf,
                    gene.name,
                    gene.n,
                    numpy.mean(muA1_post),
                    numpy.mean(muA2_post),
                    numpy.mean(muB1_post),
                    numpy.mean(muB2_post),
                    mean_logFC_A,
                    mean_logFC_B,
                    mean_delta_logFC,
                    l_delta_logFC,
                    u_delta_logFC,
                    probROPE,
                    not_HDI_overlap_bit,
                )
            )

            percentage = (100.0 * (count + 1) / N)
            progress_update(f"Running Anova Method... {percentage:5.1f}%", percentage)

            transit_tools.log(
                "analyzing %s (%1.1f%% done)" % (gene.orf, 100.0 * count / (N - 1))
            )
            count += 1

        # for HDI, maybe I should sort on abs(mean_delta_logFC); however, need to sort by prob to calculate BFDR
        probcol = -2  # probROPEs
        results.sort(key=lambda x: x[probcol])
        sortedprobs = numpy.array([x[probcol] for x in results])

        # BFDR method: Newton et al (2004). Detecting differential gene expression with a semiparametric hierarchical mixture method.  Biostatistics, 5:155-176.

        if self.inputs.signif == "BFDR":
            sortedprobs = numpy.array(sortedprobs)
            # sortedprobs.sort() # why, since already sorted?
            bfdr = numpy.cumsum(sortedprobs) / numpy.arange(1, len(sortedprobs) + 1)
            adjusted_prob = bfdr  # should be same order as sorted above by probROPE
            adjusted_label = "BFDR"

        elif self.inputs.signif == "FWER":
            fwer = stat_tools.fwer_bayes(sortedprobs)
            # fwer.sort() # should not need this if monotonic
            adjusted_prob = fwer
            adjusted_label = "FWER"

        # If not using adjustment for classification, sort correctly
        else:
            adjusted_prob = sortedprobs
            adjusted_label = "un"
            # should I stable-sort by overlap_bit?

        #            sorted_index = numpy.argsort([d[-1] for d in data])[::-1][:len(data)]
        #            adjusted_prob = [adjusted_prob[ii] for ii in sorted_index]
        #            data = [data[ii] for ii in sorted_index]

        # determine type of interactions
        extended_results = []
        for i, row in enumerate(results):
            (
                orf,
                name,
                n,
                mean_muA1_post,
                mean_muA2_post,
                mean_muB1_post,
                mean_muB2_post,
                mean_logFC_A,
                mean_logFC_B,
                mean_delta_logFC,
                l_delta_logFC,
                u_delta_logFC,
                probROPE,
                not_HDI_overlap_bit,
            ) = row

            interaction = self.classify_interaction(
                mean_delta_logFC, mean_logFC_B, mean_logFC_A
            )
            type_of_interaction = "No Interaction"
            if self.inputs.signif in "prob BFDR FWER" and adjusted_prob[i] < 0.05:
                type_of_interaction = interaction
            if self.inputs.signif == "HDI" and not_HDI_overlap_bit:
                type_of_interaction = interaction

            # switch order of last 2 vals, and append 2 more...
            new_row = tuple(
                list(row[:-2])
                + [not_HDI_overlap_bit, probROPE, adjusted_prob[i], type_of_interaction]
            )
            extended_results.append(new_row)

        return extended_results,adjusted_label # results: one row for each gene

    def print_GI(self,results,adjusted_label,condA1,condA2,condB1,condB2,sample_table): 

        self.output = sys.stdout # print to console if not output file defined
        if self.inputs.output_path != None:
            self.output = open(self.inputs.output_path, "w")
        self.output.write("%s\n" % self.identifier)

        if True: # was 'if self.wxobj:', try universal.interface=="gui" or "console" ###?
            members = sorted(
                [
                    attr
                    for attr in dir(self)
                    if not callable(getattr(self, attr)) and not attr.startswith("__")
                ]
            )
            #memberstr = ""
            #for m in members:
            #    memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write( 
                "#parameters: norm=%s, samples=%s, includeZeros=%s, output=%s\n"
                % (
                    self.inputs.normalization,
                    self.inputs.samples,
                    self.inputs.includeZeros,
                    self.output.name.encode("utf-8"), ###? add more params
                )
            )
        # originally, we wrote CLI args from console into output file, if console mode:
        #    self.output.write("#Console: python3 %s\n" % " ".join(sys.argv)) ###? I should print command

        now = str(datetime.datetime.now())
        now = now[: now.rfind(".")]
        self.output.write("#Date: " + now + "\n")
        # self.output.write("#Runtime: %s s\n" % (time.time() - start_time))

        condA1samples = [x["name"] for x in sample_table.rows if x["condition"]==condA1]
        condA2samples = [x["name"] for x in sample_table.rows if x["condition"]==condA2]
        condB1samples = [x["name"] for x in sample_table.rows if x["condition"]==condB1]
        condB2samples = [x["name"] for x in sample_table.rows if x["condition"]==condB2]

        self.output.write("#Strain A, Condition 1): %s: %s\n" % (condA1,','.join(condA1samples)))
        self.output.write("#Strain A, Condition 2): %s: %s\n" % (condA2,','.join(condA2samples)))
        self.output.write("#Strain B, Condition 1): %s: %s\n" % (condB1,','.join(condB1samples)))
        self.output.write("#Strain B, Condition 2): %s: %s\n" % (condB2,','.join(condB2samples)))
        self.output.write("#Annotation path: %s\n" % (self.inputs.annotation_path.encode("utf-8")))
        self.output.write("#ROPE=%s, method for significance=%s\n" % (self.inputs.rope, self.inputs.signif))
        if self.inputs.signif == "HDI":
            self.output.write("#Significant interactions are those genes whose delta-logFC HDI does not overlap the ROPE\n")
        elif self.inputs.signif in "prob BDFR FWER":
            self.output.write("#Significant interactions are those whose %s-adjusted probability of the delta-logFC falling within ROPE is < 0.05.\n" % (adjusted_label))

        # Write column names (redundant with self.columns)
        self.output.write(
            "#ORF\tName\tNumber of TA Sites\tMean count (Strain A Condition 1)\tMean count (Strain A Condition 2)\tMean count (Strain B Condition 1)\tMean count (Strain B Condition 2)\tMean logFC (Strain A)\tMean logFC (Strain B) \tMean delta logFC\tLower Bound delta logFC\tUpper Bound delta logFC\tIs HDI outside ROPE?\tProb. of delta-logFC being within ROPE\t%s-Adjusted Probability\tType of Interaction\n"
            % adjusted_label
        )

        # Write gene results
        for i, new_row in enumerate(results):
            (
                orf,
                name,
                n,
                mean_muA1_post,
                mean_muA2_post,
                mean_muB1_post,
                mean_muB2_post,
                mean_logFC_A,
                mean_logFC_B,
                mean_delta_logFC,
                l_delta_logFC,
                u_delta_logFC,
                not_HDI_overlap_bit,
                probROPE,
                adjusted_prob,
                type_of_interaction
            ) = new_row

            self.output.write(
                "%s\t%s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%s\t%1.8f\t%1.8f\t%s\n"
                % new_row
            )

        if self.inputs.output_path != None: self.output.close() # otherwise, it is sys.stdout
        transit_tools.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)                         
        #self.finish()
        transit_tools.log("Finished Genetic Interactions Method")

    @staticmethod
    def classify_interaction(delta_logFC, logFC_KO, logFC_WT):
        if delta_logFC < 0:
            return "Aggravating"
        elif delta_logFC >= 0 and abs(logFC_KO) < abs(logFC_WT):
            return "Alleviating"
        elif delta_logFC >= 0 and abs(logFC_KO) > abs(logFC_WT):
            return "Suppressive"
        else:
            return "N/A"


@transit_tools.ResultsFile
class File(Analysis):
    @staticmethod
    def can_load(path):
        with open(path) as in_file:
            for line in in_file:
                if line.startswith("#"):
                    if line.startswith(Analysis.identifier):
                        return True
                else:
                    return False
        return False
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Analysis.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title=self.short_desc,heading="",column_names=self.column_names,rows=self.rows).Show(),
            })
        )
        
        # 
        # get column names
        # 
        comments, headers, rows = csv.read(self.path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
        if len(comments) == 0:
            raise Exception(f'''No comments in file, and I expected the last comment to be the column names, while to load GI file "{self.path}"''')
        self.column_names = comments[-1].split("\t")
        
        # 
        # get rows
        #
        self.rows = []
        for each_row in rows:
            row = {}
            for each_column_name, each_cell in zip(self.column_names, each_row):
               row[each_column_name] = each_cell
            self.rows.append(row)
        
    
    def __str__(self):
        return f"""
            File for {self.short_name}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
    
    
Method = GUI = Analysis
Analysis() # make sure there's one instance
