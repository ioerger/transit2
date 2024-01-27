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

from pytransit.specific_tools.transit_tools import wx, basename
from pytransit.specific_tools import  gui_tools, transit_tools, console_tools, tnseq_tools, norm_tools
from pytransit.generic_tools import csv, misc
import pytransit.components.file_display as file_display
import pytransit.components.results_area as results_area
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.components.parameter_panel import panel, progress_update
from pytransit.components.spreadsheet import SpreadSheet

from pytransit.methods.pathway_enrichment import Method as PathwayEnrichment

@misc.singleton
class Method:
    name = "Genetic Interaction"
    identifier  = "GI"
    cli_name    = identifier.lower()
    menu_name   = f"{identifier} - Genetic Interaction analysis"
    description = """Genetic Interaction analysis"""
    transposons = [ "himar1" ]
    
    significance_threshold = 0.05
    
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
        "--n", # normalization flag
        "--s", # number of samples
        "--iN", # trimming TA sites at gene termini
        "--iC", 
        "--signif", # method to determine genes with significant interaction
        "--rope", # Region Of Probable Equivalence around 0
    ]

    usage_string = f"""
        Usage: {console_tools.subcommand_prefix} {cli_name} <combined_wig_file> <metadata_file> <annotation_file> <conditionA1> <conditionB1> <conditionA2> <conditionB2> <output_file> [optional arguments]
            GI performs a comparison among 2x2=4 groups of datasets, e.g. strains A and B assessed in conditions 1 and 2 (e.g. control vs treatment).
            It looks for interactions where the response to the treatment (i.e. effect on insertion counts) depends on the strain (output variable: delta_LFC).
            Provide replicates in each group as a comma-separated list of wig files.
            HDI is highest density interval for posterior distribution of delta_LFC, which is like a confidence interval on difference of slopes.
            Genes are sorted by probability of HDI overlapping with ROPE. (genes with the highest abs(mean_delta_logFC) are near the top, approximately)
            Significant genes are indicated by 'Type of Interaction' column (No Interaction, Aggravating, Alleviating, Suppressive).
                By default, hits are defined as "Is HDI outside of ROPE?"=TRUE (i.e. non-overlap of delta_LFC posterior distritbuion with Region of Probably Equivalence around 0)
                Alternative methods for significance: use --signif key with prob, BFDR, or FWER. These affect 'Type of Interaction' (i.e. which genes are labeled 'No Interaction')

        Optional Arguments:
            --n <string>     :=  Normalization method. Default: --n TTR
            --s <integer>    :=  Number of samples. Default: --s 10000
            --iN <float>     :=  Ignore TAs occurring at given percentage (as integer) of the N terminus. Default: --iN 0
            --iC <float>     :=  Ignore TAs occurring at given percentage (as integer) of the C terminus. Default: --iC 0
            --rope <float>   :=  Region of Practical Equivalence. Area around 0 (i.e. 0 +/- ROPE) that is NOT of interest. Can be thought of similar to the area of the null-hypothesis. Default: --rope 0.5
            --signif HDI     :=  (default) Significant if HDI does not overlap ROPE; if HDI overlaps ROPE, 'Type of Interaction' is set to 'No Interaction'
            --signif prob    :=  Optionally, significant hits are re-defined based on probability (degree) of overlap of HDI with ROPE, prob<{significance_threshold} (no adjustment)
            --signif BFDR    :=  Apply "Bayesian" FDR correction (see doc) to adjust HDI-ROPE overlap probabilities so that significant hits are re-defined as BFDR<{significance_threshold}
            --signif FWER    :=  Apply "Bayesian" FWER correction (see doc) to adjust HDI-ROPE overlap probabilities so that significant hits are re-defined as FWER<{significance_threshold}
    """.replace("\n        ", "\n")
    
    @gui.add_menu("Method", "Himar1", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    #@gui.add_menu("Method", "Tn5", menu_name)
    #def on_menu_click(event):
    #    Method.define_panel(event)

    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            set_instructions(
                title_text= self.name,
                sub_text = self.description,
                method_specific_instructions="""
                    GI performs a comparison among 2x2=4 groups of datasets, e.g. strains A and B assessed in conditions 1 and 2 (e.g. control vs treatment). It looks for interactions where the response to the treatment (i.e. effect on insertion counts) depends on the strain (output variable: delta_LFC). Provide replicates in each group as a comma-separated list of wig files.
                    
                    HDI is highest density interval for posterior distribution of delta_LFC, which is like a confidence interval on difference of slopes. Genes are sorted by probability of HDI overlapping with ROPE. (genes with the highest abs(mean_delta_logFC) are near the top, approximately). Significant genes are indicated by 'Type of Interaction' column (No Interaction, Aggravating, Alleviating, Suppressive).
                    
                    By default, hits are defined as "Is HDI outside of ROPE?"=TRUE (i.e. non-overlap of delta_LFC posterior distritbuion with Region of Probably Equivalence around 0)

                    1.  Add an annotation file for the organism corresponding to the desired datasets

                    2.  Adding datasets grown under condition A
                        a. Using the dropdown for Condition A1, select your control dataset for condition A
                        b. Using the dropdown for Condition A2, select your experimental dataset for condition A

                    3.  Adding datasets grown under condition B
                        a. Using the dropdown for Condition B1, select your control dataset for condition B
                        b. Using the dropdown for Condition B2, select your experimental dataset for condition B

                    4. [Optional] Select the remaining parameters

                    5. Click Run
                """.replace("\n                    ","\n"),
            )

            # only need Norm selection and Run button        
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
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


    @staticmethod
    def from_gui(frame):
        # 
        # get wig files
        # 
        wig_group = gui.combined_wigs[-1] # assume there is only 1 (should check that it has beed defined)
        Method.inputs.combined_wig_path = wig_group.main_path # see components/sample_area.py
        Method.inputs.metadata_path = gui.combined_wigs[-1].metadata_path # assume all samples are in the same metadata file

        # 
        # get annotation
        # 
        Method.inputs.annotation_path = gui.annotation_path
        
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
            default_file_name=f"{Method.cli_name}_output.tsv",
            output_extensions=transit_tools.result_output_extensions,
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

        combined_wig_path = args[0]
        metadata_path     = args[1]
        annotation_path   = args[2]
        condA1            = args[3]
        condA2            = args[4]
        condB1            = args[5]
        condB2            = args[6]
        output_path       = args[7]

        normalization = kwargs.get("n", "TTR")
        samples       = int(kwargs.get("s", 10000))
        rope          = float(kwargs.get("rope", 0.5))  # fixed! changed int to float
        signif        = kwargs.get("signif", "HDI")
        n_terminus    = float(kwargs.get("iN", 0.00))
        c_terminus    = float(kwargs.get("iC", 0.00))

        # save all the data
        Method.inputs.update(dict(
          combined_wig_path=combined_wig_path,
          annotation_path=annotation_path,
          metadata_path=metadata_path,
          condA1=condA1,
          condA2=condA2,
          condB1=condB1,
          condB2=condB2,
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
        with gui_tools.nice_error_log:
            logging.log("Starting Genetic Interaction analysis")
            self.start_time = time.time()

            ##########################
            # get data

            logging.log("Getting Data")
            sites, data, filenames_in_comb_wig = tnseq_tools.CombinedWigData.load(self.inputs.combined_wig_path)
            logging.log(f"Normalizing using: {self.inputs.normalization}")
            data, factors = norm_tools.normalize_data(data, self.inputs.normalization)

            # is it better to read the metadata directly, rather than pulling from samples_table, to accommodate console mode?
            metadata = tnseq_tools.CombinedWigMetadata(self.inputs.metadata_path)


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
            (results,adjusted_label) = self.calc_gi(dataA1,dataA2,dataB1,dataB2,sites)

            ######################
            # write output
            # 
            # note: first comment line is filetype, last comment line is column headers

            logging.log(f"Adding File: {self.inputs.output_path}")
            
            # will open and close file, or print to console
            self.print_gi_results(results,adjusted_label,condA1,condA2,condB1,condB2,metadata)
            

            aggra = len(list(filter(lambda x: x[-1]=="Aggravating", results)))
            allev = len(list(filter(lambda x: x[-1]=="Alleviating", results)))
            suppr = len(list(filter(lambda x: x[-1]=="Suppressive", results)))
            logging.log("Summary of genetic interactions: aggravating=%s, alleviating=%s, suppressive=%s" % (aggra,allev,suppr))
            logging.log("Time: %0.1fs\n" % (time.time() - self.start_time))
            logging.log("Finished Genetic Interaction analysis")
            results_area.add(self.inputs.output_path)


    def calc_gi(self, dataA1, dataA2, dataB1, dataB2, position): # position is vector of TAsite coords
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

                alpha = Method.significance_threshold # TODO: not sure if these two are equivlent by coincidence or by nature

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
            progress_update(f"Running GI Method... {percentage:5.1f}%", percentage)

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
            if self.inputs.signif in "prob BFDR FWER" and adjusted_prob[i] < Method.significance_threshold:
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

    def print_gi_results(self,results,adjusted_label,condA1,condA2,condB1,condB2,metadata): 
        # 
        # computation
        # 
        condA1samples = [x["Id"] for x in metadata.rows if x["Condition"]==condA1]
        condA2samples = [x["Id"] for x in metadata.rows if x["Condition"]==condA2]
        condB1samples = [x["Id"] for x in metadata.rows if x["Condition"]==condB1]
        condB2samples = [x["Id"] for x in metadata.rows if x["Condition"]==condB2]
        
        if self.inputs.signif == "HDI":
            significance_note = "Significant interactions are those genes whose delta-logFC HDI does not overlap the ROPE"
        elif self.inputs.signif in "prob BDFR FWER":
            significance_note = f"Significant interactions are those whose {adjusted_label}-adjusted probability of the delta-logFC falling within ROPE is < {Method.significance_threshold}."
        else:
            significance_note = None
        
        aggra = len(list(filter(lambda x: x[-1]=="Aggravating", results)))
        allev = len(list(filter(lambda x: x[-1]=="Alleviating", results)))
        suppr = len(list(filter(lambda x: x[-1]=="Suppressive", results)))
        
        annotation_data = tnseq_tools.AnnotationFile(path=self.inputs.annotation_path)
        
        # 
        # format into rows
        # 
        rows = []        
        for orf, name, n, mean_mu_a1_post, mean_mu_a2_post, mean_mu_b1_post, mean_mu_b2_post, mean_log_fc_a, mean_log_fc_b, mean_delta_log_fc, l_delta_log_fc, u_delta_log_fc, not_hdi_overlap_bit, prob_rope, adjusted_prob, type_of_interaction in results:
            rows.append([
                orf,
                name,
                annotation_data.gene_description(orf_id=orf),
                "%d" % n,
                "%1.2f" % mean_mu_a1_post,
                "%1.2f" % mean_mu_a2_post,
                "%1.2f" % mean_mu_b1_post,
                "%1.2f" % mean_mu_b2_post,
                "%1.2f" % mean_log_fc_a,
                "%1.2f" % mean_log_fc_b,
                "%1.2f" % mean_delta_log_fc,
                "%1.2f" % l_delta_log_fc,
                "%1.2f" % u_delta_log_fc,
                str(not_hdi_overlap_bit),
                "%1.8f" % prob_rope,
                "%1.8f" % adjusted_prob,
                str(type_of_interaction),
            ])
        
        # 
        # actual write
        # 
        transit_tools.write_result(
            path=self.inputs.output_path,
            file_kind=Method.identifier,
            rows=rows,
            column_names=[
                'ORF',
                'Gene',
                'Annotation',
                'TA Sites',
                'A1 Mean Count',
                'A2 Mean Count',
                'B1 Mean Count',
                'B2 Mean Count',
                'Log 2 FC Strain A',
                'Log 2 FC Strain B',
                'Delta Log 2 FC',
                'Lower Bound Delta Log 2 FC',
                'Upper Bound Delta Log 2 FC',
                'Is HDI Outside ROPE',
                'Probability Of Delta Log 2 FC Being Within ROPE',
                str(adjusted_label)+ ' Adj P Value',
                'Type Of Interaction',
            ],
            extra_info=dict(
                calculation_time=f"{time.time() - self.start_time:0.1f}seconds",
                analysis_type=Method.identifier,
                conditions=", ".join([condA1,condA1,condB1,condB2]),
                files=dict(
                    combined_wig=Method.inputs.combined_wig_path,
                    annotation_path=Method.inputs.annotation_path,
                ),
                parameters= dict(
                    Normalization_Of_Counts= self.inputs.normalization,
                    Number_Of_Samples_For_Monte_Carlo = str(self.inputs.samples),
                    Trimming_Of_TA_Sites = dict(
                        N_Terminus = str(self.inputs.n_terminus),
                        C_Terminus = str(self.inputs.c_terminus),
                    ),
                    ROPE = str(self.inputs.rope)+" (region of probable equivalence around 0)",
                    Method_For_Determining_Significance = str(self.inputs.signif),
                    Annotation_Path = self.inputs.annotation_path,
                    Output_Path = self.inputs.output_path,
                ),
                Samples_In_4_Conditions= dict(
                    Strain_A_Condition_1= str(condA1)+" = " +','.join([str(x) for x in condA1samples]),
                    Strain_A_Condition_2= str(condA2)+" = " +','.join([str(x) for x in condA2samples]),
                    Strain_B_Condition_1= str(condB1)+" = " + ','.join([str(x) for x in condB1samples]),
                    Strain_B_Condition_2= str(condB2)+" = " + ','.join([str(x) for x in condB2samples]),
                ),
                Significance_Note = significance_note,
                Summary_Of_Genetic_Interactions = dict(
                    Aggravating = str(aggra),
                    Alleviating = str(allev),
                    Suppressive = str(suppr),
                ),
            ),
        )
        
        logging.log("Adding File: %s" % (self.inputs.output_path))                       
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
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.description,
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows
                    ).Show(),
                    "Pathway Enrichment": lambda *args: PathwayEnrichment.call_from_results_panel(path),
            })
        )
        
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        summary = self.extra_data["Summary_Of_Genetic_Interactions"]
        summary_str = [str(summary[key])+" "+str(key) for key in sorted(summary.keys())] 
        self.values_for_result_table.update({"summary": "; ".join(summary_str) })
        #self.values_for_result_table.update(self.extra_data.get("summary_info", {}))

        parameters = self.extra_data.get("parameters",{})
        parameters_str = [str(key)+" : "+str(parameters[key]) for key in ["Number_Of_Samples_For_Monte_Carlo","Normalization_Of_Counts","ROPE","Method_For_Determining_Significance"]]
        self.values_for_result_table.update({"parameters": "; ".join(parameters_str) })
        
       
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
    
    
