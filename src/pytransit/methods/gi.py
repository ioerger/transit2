from pytransit.components.parameter_panel import panel, progress_update, set_instructions
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
from pytransit.generic_tools.lazy_dict import LazyDict

from pytransit.specific_tools.transit_tools import wx, basename, HAS_R, FloatVector, DataFrame, StrVector
from pytransit.specific_tools import logging, gui_tools, transit_tools, console_tools, tnseq_tools, norm_tools
from pytransit.generic_tools import csv, misc
import pytransit.components.file_display as file_display
import pytransit.components.results_area as results_area
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components.parameter_panel import panel, progress_update
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    name = "Genetic Interaction"
    identifier  = "GI"
    cli_name    = identifier.lower()
    menu_name   = f"{identifier} - Genetic Interaction analysis"
    description = """Genetic Interaction analysis"""
    transposons = [ "himar1" ]
    
    inputs = LazyDict(
          combined_wig=None,
          metadata_path=None,
          condA1=None,
          condA2=None,
          condB1=None,
          condB2=None,

          annotation_path=None,
          output_path=None,

          normalization=None,
          samples=None,
          rope=None,
          signif=None,

          n_terminus=None,
          c_terminus=None,
    )
    
    valid_cli_flags = [
        "-n", # normalization flag
        "-s", # number of samples
        "-iN", # trimming TA sites at gene termini
        "-iC", 
        "-signif", # method to determine genes with significant interaction
        "--rope", # Region Of Probable Equivalence around 0
    ]

    usage_string = f"""usage: {console_tools.subcommand_prefix} gi <combined_wig> <samples_metadata> <conditionA1> <conditionB1> <conditionA2> <conditionB2> <prot_table> <output_file> [optional arguments]
        GI performs a comparison among 2x2=4 groups of datasets, e.g. strains A and B assessed in conditions 1 and 2 (e.g. control vs treatment).
        It looks for interactions where the response to the treatment (i.e. effect on insertion counts) depends on the strain (output variable: delta_LFC).
        Provide replicates in each group as a comma-separated list of wig files.
        HDI is highest density interval for posterior distribution of delta_LFC, which is like a confidence interval on difference of slopes.
        Genes are sorted by probability of HDI overlapping with ROPE. (genes with the highest abs(mean_delta_logFC) are near the top, approximately)
        Significant genes are indicated by 'Type of Interaction' column (No Interaction, Aggravating, Alleviating, Suppressive).
            By default, hits are defined as "Is HDI outside of ROPE?"=TRUE (i.e. non-overlap of delta_LFC posterior distritbuion with Region of Probably Equivalence around 0)
            Alternative methods for significance: use -signif flag with prob, BFDR, or FWER. These affect 'Type of Interaction' (i.e. which genes are labeled 'No Interaction')

        Optional Arguments:
        -n <string>     :=  Normalization method. Default: -n TTR
        -s <integer>    :=  Number of samples. Default: -s 10000
        -iN <float>     :=  Ignore TAs occuring at given percentage (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring at given percentage (as integer) of the C terminus. Default: -iC 0
        --rope <float>  :=  Region of Practical Equivalence. Area around 0 (i.e. 0 +/- ROPE) that is NOT of interest. Can be thought of similar to the area of the null-hypothesis. Default: --rope 0.5
        -signif HDI     :=  (default) Significant if HDI does not overlap ROPE; if HDI overlaps ROPE, 'Type of Interaction' is set to 'No Interaction'
        -signif prob    :=  Optionally, significant hits are re-defined based on probability (degree) of overlap of HDI with ROPE, prob<0.05 (no adjustment)
        -signif BFDR    :=  Apply "Bayesian" FDR correction (see doc) to adjust HDI-ROPE overlap probabilities so that significant hits are re-defined as BFDR<0.05
        -signif FWER    :=  Apply "Bayesian" FWER correction (see doc) to adjust HDI-ROPE overlap probabilities so that significant hits are re-defined as FWER<0.05
    """
    
    @gui.add_menu("Method", "himar1", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    @gui.add_menu("Method", "tn5", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)

    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            set_instructions(
                method_short_text= self.name,
                method_long_text = self.description,
                method_descr="""
                    GI performs a comparison among 2x2=4 groups of datasets, e.g. strains A and B assessed in conditions 1 and 2 (e.g. control vs treatment).
                    It looks for interactions where the response to the treatment (i.e. effect on insertion counts) depends on the strain (output variable: delta_LFC).
                    Provide replicates in each group as a comma-separated list of wig files.
                    HDI is highest density interval for posterior distribution of delta_LFC, which is like a confidence interval on difference of slopes.
                    Genes are sorted by probability of HDI overlapping with ROPE. (genes with the highest abs(mean_delta_logFC) are near the top, approximately)
                    Significant genes are indicated by 'Type of Interaction' column (No Interaction, Aggravating, Alleviating, Suppressive).
                    By default, hits are defined as "Is HDI outside of ROPE?"=TRUE (i.e. non-overlap of delta_LFC posterior distritbuion with Region of Probably Equivalence around 0)
                    Alternative methods for significance: use -signif flag with prob, BFDR, or FWER. These affect 'Type of Interaction' (i.e. which genes are labeled 'No Interaction')
                """.replace("\n            ","\n"),
                method_specific_instructions="""
                    FIXME
                """.replace("\n            ","\n")
            )

            # only need Norm selection and Run button        
            self.value_getters = LazyDict()
            self.value_getters.condA1        = panel_helpers.create_condition_input(panel,main_sizer, label_text="Condition A1:", tooltip_text="indicate condition representing 'strain A' in 'condition 1'")
            self.value_getters.condB1        = panel_helpers.create_condition_input(panel,main_sizer, label_text="Condition B1:", tooltip_text="indicate condition representing 'strain B' in 'condition 1'")
            self.value_getters.condA2        = panel_helpers.create_condition_input(panel,main_sizer, label_text="Condition A2:", tooltip_text="indicate condition representing 'strain A' in 'condition 2'")
            self.value_getters.condB2        = panel_helpers.create_condition_input(panel,main_sizer, label_text="Condition B2:", tooltip_text="indicate condition representing 'strain B' in 'condition 2'")
            self.value_getters.normalization = panel_helpers.create_normalization_input(panel, main_sizer) # TTR is default
            self.value_getters.n_terminus    = panel_helpers.create_n_terminus_input(panel, main_sizer)
            self.value_getters.c_terminus    = panel_helpers.create_c_terminus_input(panel, main_sizer)
            self.value_getters.samples       = panel_helpers.create_int_getter(  panel, main_sizer, label_text="Number of samples",default_value=10000, tooltip_text="random trials in Monte Carlo simulation")
            self.value_getters.rope          = panel_helpers.create_float_getter(panel, main_sizer, label_text="ROPE"             ,default_value=0.5,   tooltip_text="Region of probable equivalence around 0")
            self.value_getters.signif        = panel_helpers.create_significance_choice_box(panel,main_sizer) # default is HDI
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)


    @staticmethod
    def from_gui(frame):
        # 
        # get wig files
        # 
        wig_group = gui.combined_wigs[0] # assume there is only 1 (should check that it has beed defined)
        Method.inputs.combined_wig = wig_group.main_path # see components/sample_area.py
        Method.inputs.metadata_path = gui.combined_wigs[0].metadata_path # assume all samples are in the same metadata file

        # 
        # get annotation
        # 
        Method.inputs.annotation_path = gui.annotation_path
        transit_tools.validate_annotation(Method.inputs.annotation_path)
        
        # 
        # setup custom inputs
        # 
        for each_key, each_getter in Method.value_getters.items():
            try:
                Method.inputs[each_key] = each_getter()
            except Exception as error:
                raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')

        # 
        # save result files
        # 
        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name="GI_test.dat",
            output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
        )
        if not Method.inputs.output_path:
            return None

        return Method

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=8)

        combined_wig = args[0]
        metadata_path = args[1]
        condA1 = args[2]
        condA2 = args[3]
        condB1 = args[4]
        condB2 = args[5]

        annotation_path = args[6]
        output_path = args[7]

        normalization = kwargs.get("n", "TTR")
        samples = int(kwargs.get("s", 10000))
        rope = float(kwargs.get("-rope", 0.5))  # fixed! changed int to float
        signif = kwargs.get("signif", "HDI")

        n_terminus = float(kwargs.get("iN", 0.00))
        c_terminus = float(kwargs.get("iC", 0.00))

        # save all the data
        Method.inputs.update(dict(
          combined_wig=combined_wig,
          metadata_path=metadata_path,
          condA1=condA1,
          condA2=condA2,
          condB1=condB1,
          condB2=condB2,

          annotation_path=annotation_path,
          output_path=output_path,

          normalization=normalization,
          samples=samples,
          rope=rope,
          signif=signif,

          n_terminus=n_terminus,
          c_terminus=c_terminus
        ))
        
        Method.Run()
        
    def Run(self):
        logging.log("Starting Genetic Interaction analysis")
        start_time = time.time()

        ##########################
        # get data

        logging.log("Getting Data")
        sites, data, filenames_in_comb_wig = tnseq_tools.read_combined_wig(self.inputs.combined_wig)
        logging.log(f"Normalizing using: {self.inputs.normalization}")
        data, factors = norm_tools.normalize_data(data, self.inputs.normalization)
            
        # # Do LOESS correction if specified
        # if self.LOESS:
        #   logging.log("Performing LOESS Correction")
        #   for j in range(K):
        #     data[j] = stat_tools.loess_correction(position, data[j])

        # is it better to read the metadata directly, rather than pulling from samples_table, to accommodate console mode?
        metadata = tnseq_tools.CombinedWigMetadata(self.inputs.metadata_path)
        #for sample in metadata.rows:
        #  print("%s\t%s" % (sample["Id"],sample["Condition"]))


        ##########################
        # process data
        # 
        # note: self.inputs.condA1 = self.value_getters.condA1()
        logging.log("processing data")
        # get 4 lists of indexes into data (extract columns for 4 conds in comwig)

        indexes = {}
        for i,row in enumerate(metadata.rows): 
            cond = row["Condition"] # "condition" for samples_table
            if cond not in indexes: indexes[cond] = []
            indexes[cond].append(i) 
        condA1 = self.inputs.condA1
        condA2 = self.inputs.condA2
        condB1 = self.inputs.condB1
        condB2 = self.inputs.condB2

        if condA1 not in indexes or len(indexes[condA1])==0: logging.error("no samples found for condition %s" % condA1)
        if condA2 not in indexes or len(indexes[condA2])==0: logging.error("no samples found for condition %s" % condA2)
        if condB1 not in indexes or len(indexes[condB1])==0: logging.error("no samples found for condition %s" % condB1)
        if condB2 not in indexes or len(indexes[condB2])==0: logging.error("no samples found for condition %s" % condB2)

        logging.log("condA1=%s, samples=%s" % (condA1,','.join([str(x["Id"]) for x in metadata.rows if x["Condition"]==condA1])))
        logging.log("condA2=%s, samples=%s" % (condA1,','.join([str(x["Id"]) for x in metadata.rows if x["Condition"]==condA2])))
        logging.log("condB1=%s, samples=%s" % (condA1,','.join([str(x["Id"]) for x in metadata.rows if x["Condition"]==condB1])))
        logging.log("condB2=%s, samples=%s" % (condA1,','.join([str(x["Id"]) for x in metadata.rows if x["Condition"]==condB2])))

        dataA1 = data[indexes[condA1]] # select datasets (rows)
        dataA2 = data[indexes[condA2]]
        dataB1 = data[indexes[condB1]]
        dataB2 = data[indexes[condB2]]

        # results: 1 row for each gene; adjusted_label - just a string that gets printed in the header
        (results,adjusted_label) = self.calc_GI(dataA1,dataA2,dataB1,dataB2,sites)

        ######################
        # write output
        # 
        # note: first comment line is filetype, last comment line is column headers

        logging.log(f"Adding File: {self.inputs.output_path}")
        results_area.add(self.inputs.output_path)

        # will open and close file, or print to console
        self.print_GI_results(results,adjusted_label,condA1,condA2,condB1,condB2,metadata)

        aggra = len(list(filter(lambda x: x[-1]=="Aggravating", results)))
        allev = len(list(filter(lambda x: x[-1]=="Alleviating", results)))
        suppr = len(list(filter(lambda x: x[-1]=="Suppressive", results)))
        logging.log("Summary of genetic interactions: aggravating=%s, alleviating=%s, suppressive=%s" % (aggra,allev,suppr))
        logging.log("Time: %0.1fs\n" % (time.time() - start_time))
        logging.log("Finished Genetic Interaction analysis")


    def calc_GI(self,dataA1,dataA2,dataB1,dataB2,position): # position is vector of TAsite coords
        from pytransit.specific_tools import stat_tools
        
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

            logging.log(
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

        return extended_results, adjusted_label # results: one row for each gene

    def print_GI_results(self,results,adjusted_label,condA1,condA2,condB1,condB2,metadata): 

        self.output = sys.stdout # print to console if not output file defined
        if self.inputs.output_path != None:
            self.output = open(self.inputs.output_path, "w")
        self.output.write("#%s\n" % self.identifier)

        if not gui.is_active:
            self.output.write(f"#Console: {console_tools.full_commandline_command}")

        now = str(datetime.datetime.now())
        now = now[: now.rfind(".")]
        self.output.write("#Date: " + now + "\n")
        # self.output.write("#Runtime: %s s\n" % (time.time() - start_time))

        self.output.write("#Parameters:\n")
        self.output.write("# normalization of counts=%s\n" % (self.inputs.normalization))
        self.output.write("# num samples for Monte Carlo=%s\n" % (self.inputs.samples))
        self.output.write("# trimming of TA sites: Nterm=%s%%, Cterm=%s%%\n" % (self.inputs.n_terminus,self.inputs.c_terminus))
        self.output.write("# ROPE=%s (region of probable equivalence around 0)\n" % (self.inputs.rope))
        self.output.write("# method for determining significance=%s\n" % (self.inputs.signif))
        self.output.write("# annotation path: %s\n" % (self.inputs.annotation_path.encode("utf-8")))
        self.output.write("# output path: %s\n" % (self.inputs.output_path))

        condA1samples = [x["Id"] for x in metadata.rows if x["Condition"]==condA1]
        condA2samples = [x["Id"] for x in metadata.rows if x["Condition"]==condA2]
        condB1samples = [x["Id"] for x in metadata.rows if x["Condition"]==condB1]
        condB2samples = [x["Id"] for x in metadata.rows if x["Condition"]==condB2]

        self.output.write("#Samples in 4 conditions:\n")
        self.output.write("# A1 (Strain A, Condition 1): %s: %s\n" % (condA1,','.join([str(x) for x in condA1samples])))
        self.output.write("# A2 (Strain A, Condition 2): %s: %s\n" % (condA2,','.join([str(x) for x in condA2samples])))
        self.output.write("# B1 (Strain B, Condition 1): %s: %s\n" % (condB1,','.join([str(x) for x in condB1samples])))
        self.output.write("# B2 (Strain B, Condition 2): %s: %s\n" % (condB2,','.join([str(x) for x in condB2samples])))
        if self.inputs.signif == "HDI":
            self.output.write("#Significant interactions are those genes whose delta-logFC HDI does not overlap the ROPE\n")
        elif self.inputs.signif in "prob BDFR FWER":
            self.output.write("#Significant interactions are those whose %s-adjusted probability of the delta-logFC falling within ROPE is < 0.05.\n" % (adjusted_label))

        aggra = len(list(filter(lambda x: x[-1]=="Aggravating", results)))
        allev = len(list(filter(lambda x: x[-1]=="Alleviating", results)))
        suppr = len(list(filter(lambda x: x[-1]=="Suppressive", results)))
        self.output.write("#Summary of genetic interactions: aggravating=%s, alleviating=%s, suppressive=%s\n" % (aggra,allev,suppr))

        # Write column names (redundant with self.columns)
        column_names = [
            'ORF',
            'Gene',
            'Annotation',
            'TA Sites',
            'A1 Mean Count',
            'A2 Mean Count',
            'B1 Mean Count',
            'B2 Mean Count',
            'Log FC Strain A',
            'Log FC Strain B',
            'Delta Log FC',
            'Lower Bound Delta Log FC',
            'Upper Bound Delta Log FC',
            'Is HDI Outside ROPE',
            'Probability Of Delta Log FC Being Within ROPE',
            f'{adjusted_label} Adj P Value',
            'Type Of Interaction'
        ]
        self.output.write(
            "#"+"\t".join(column_names)+"\n"
        )
    
        annot = {}
        for line in open(self.inputs.annotation_path):
          w = line.rstrip().split('\t')
          annot[w[8]] = w[0]

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
                type_of_interaction,
            ) = new_row

            orf = new_row[0]
            descr = annot[orf]
            new_row = tuple(list(new_row[:2])+[descr]+list(new_row[2:]))

            self.output.write(
                "%s\t%s\t%s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%s\t%1.8f\t%1.8f\t%s\n"
                % new_row
            )

        if self.inputs.output_path != None: self.output.close() # otherwise, it is sys.stdout
        logging.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)                         
        logging.log("Finished Genetic Interactions Method")

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
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title=Method.description,heading="",column_names=self.column_names,rows=self.rows).Show(),
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
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
    
    
