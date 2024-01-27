import sys
import os
import time
import ntpath
import math
import random
import datetime
import itertools
import statistics
import heapq
import collections
import numpy

from pytransit.generic_tools.lazy_dict import LazyDict

from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.components.parameter_panel import progress_update, set_instructions
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.specific_tools import gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools import csv, misc, informative_iterator
from pytransit.generic_tools.misc import cache
import pytransit.components.results_area as results_area
from pytransit.components import parameter_panel

@cache() # cache really only helps for testing, but it helps a lot for testing
def slow_computation(gene_one_hot_encoded, ttn_vectors, Y):
    import pandas
    import statsmodels.stats.multitest
    import statsmodels.api as sm
    X1 = pandas.concat([gene_one_hot_encoded, ttn_vectors], axis=1)
    #X1 = sm.add_constant(X1)
    return X1, sm.OLS(Y.astype(float), X1.astype(float)).fit()

@misc.singleton
class Method:
    name = "TTN Fitness"
    identifier  = name.replace(" ", "")
    cli_name    = identifier.lower()
    menu_name   = f"{name} - Analyze fitness effect of (non-essential) genes"
    description = """Analyze fitness effect of (non-essential) genes using a predictive model that corrects for the bias in Himar1 insertion probabilities based on nucleotides around each TA site"""
    transposons = [ "himar1" ] # definitely Himar1 only
    
    inputs = LazyDict(
        combined_wig = None,
        metadata = None,
        condition = None, # all reps will be combined; later, allow user to select individual wigs files
        wig_files = None,
        annotation_path = None,
        genome_path = None,
        gumbel_results_path = None,
        genes_output_path = None,
        sites_output_path = None,
        normalization = "TTR",
    )


    usage_string = f"""
        Usage: {console_tools.subcommand_prefix} {cli_name} <combined_wig_file> <metadata_file> <annotation_file> <condition to analyze> <genome .fna> <gumbel_results_file> <genes output_file> <sites output_file>
    """
    
    @gui.add_menu("Method", "Himar1", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            set_instructions(
                title_text= self.name,
                sub_text="",
                method_specific_instructions="""
                TTN-Fitness provides a method for estimating the fitness of genes in a single condition, while correcting for biases in Himar1 insertion preferences at TA sites based on surrounding nucleotides. The frequency of insertions depends on nucleotides surrounding TA sites. This model captures that effect.

                Typically with individual TnSeq datasets, Gumbel and HMM are the methods used for evaluating essentiality. Gumbel distinguishes between ES (essential) from NE (non-essential). HMM adds the GD (growth-defect; suppressed counts; mutant has reduced fitness) and GA (growth advantage; inflated counts; mutant has selective advantage) categories. Quantifying the magnitude of the fitness defect is risky because the counts at individual TA sites can be noisy. Sometimes the counts at a TA site in a gene can span a wide range of very low to very high counts. The TTN-Fitness gives a more fine-grained analysis of the degree of fitness effect by taking into account the insertion preferences of the Himar1 transposon.

                1. Add an annotation file for the organism corresponding to the desired datasets

                2. Of the Conditions in the Conditions pane, select one using the 'Condition to Analyze' dropdown

                3. Use the 'Gumbel Results File' button to upload a file of the gumbel essentiallity calls for the genes using the condition selected in the "Condition to Analyze" dropdown

                4. Use the 'Load Genome Sequence File' button to use the file dialog to upload a .fna file

                5. [Optional] Select the default name of your output files. This will be generated in your local directory

                6. Click Run
                """.replace("\n                    ","\n"),
            )
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            self.value_getters = LazyDict()

            self.value_getters.condition = panel_helpers.create_condition_input(panel, main_sizer)
            self.value_getters.gumbel_results_path = panel_helpers.create_file_input(panel, main_sizer,
                button_label="Gumbel results file",
                default_file_name="glycerol_gumbel.out",
                allowed_extensions="All files (*.*)|*.*", 
                popup_title="Choose Gumbel results file", 
                tooltip_text="Must run Gumbel first to determine which genes are essential. Note: TTN-fitness estimates fitness of NON-essential genes."
            )
            self.value_getters.genome_path = panel_helpers.create_file_input(panel, main_sizer,
                popup_title="Choose genome sequence file",
                button_label="Load genome sequence file",
                default_file_name="H37Rv.fna",
                allowed_extensions="Fasta files (*.fa;*.fna;*.fasta))|*.fa;*.fna;*.fasta",
                tooltip_text="Genome sequence file (.fna) must match annotation file (.prot_table)",
            )
            self.value_getters.output_basename = panel_helpers.create_text_box_getter(panel, main_sizer,
                label_text="Basename for output files",
                default_value="ttnfitness.test",
                tooltip_text="If X is basename, then X_genes.tsv and X_sites.tsv will be generated as output files."
            )
            self.value_getters.normalization = panel_helpers.create_normalization_input(panel, main_sizer, default=self.inputs.normalization) # TTR 
            
    
    @staticmethod
    def from_gui(frame):
        with gui_tools.nice_error_log:
            combined_wig = gui.combined_wigs[-1]
            Method.inputs.combined_wig = combined_wig.main_path
            # assume all samples are in the same metadata file
            Method.inputs.metadata_path = gui.combined_wigs[-1].metadata_path 
            Method.inputs.annotation_path = gui.annotation_path

            # 
            # call all GUI getters, puts results into respective Method.inputs key-value
            # 
            for each_key, each_getter in Method.value_getters.items():
                try:
                    Method.inputs[each_key] = each_getter()
                except Exception as error:
                    raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
            
            Method.inputs.genes_output_path = "%s.genes.tsv" % (Method.inputs.output_basename)
            Method.inputs.sites_output_path = "%s.sites.tsv" % (Method.inputs.output_basename)

            return Method

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs): # clean_args() was already called in pytransit/__main__.py
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=8)
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        Method.inputs.update(dict(
            combined_wig    = args[0],
            metadata_path   = args[1],
            annotation_path = args[2],
            condition       = args[3],
            genome_path = args[4],
            gumbel_results_path = args[5],
            genes_output_path = args[6],
            sites_output_path = args[7],
        ))
            
        Method.Run()
        
    def Run(self):
        with gui_tools.nice_error_log:
            logging.log("Starting tnseq_stats analysis")
            self.start_time = time.time()

            #######################
            # get data

            if self.inputs.combined_wig!=None:  # assume metadata and condition are defined too
                logging.log("Getting Data from %s" % self.inputs.combined_wig)
                position, data, filenames_in_comb_wig = tnseq_tools.CombinedWigData.load(self.inputs.combined_wig)

                metadata = tnseq_tools.CombinedWigMetadata(self.inputs.metadata_path)
                indexes,ids = [],[]
                for i,f in enumerate(filenames_in_comb_wig): 
                  cond = metadata.conditions_by_wig_fingerprint.get(f, "FLAG-UNMAPPED-CONDITION-IN-WIG")
                  id = metadata.id_for(f)
                  if cond==Method.inputs.condition:
                    indexes.append(i)
                    ids.append(id)
                logging.log("selected samples for ttnfitness (cond=%s): %s" % (Method.inputs.condition,','.join(ids)))
                data = data[indexes]
                # now, select the columns in data corresponding to samples that are replicates of desired condition...

            elif self.inputs.wig_files!=None:
                logging.log("Getting Data")
                (data, position) = transit_tools.get_validated_data(self.inputs.wig_files)

            else:
                logging.error("error: must provide either combined_wig or list of wig files")
                
            (K, N) = data.shape 

            # normalize the counts
            if self.inputs.normalization and self.inputs.normalization != "nonorm":
              logging.log("Normalizing using: %s" % self.inputs.normalization)
              (data, factors) = norm_tools.normalize_data( data, self.inputs.normalization, self.inputs.wig_files, self.inputs.annotation_path )
                
            # read-in genes from annotation
            G = tnseq_tools.Genes(
                self.inputs.wig_files,
                self.inputs.annotation_path,
                data=data,
                position=position,
            )
            N = len(G)
    
            logging.log("Getting Genome")
            genome = ""
            n = 0
            with open(self.inputs.genome_path) as file:
                for line in file:
                    if n == 0:
                        n = 1  # skip first
                    else:
                        genome += line[:-1]

            # could also read-in gumbel_results_file as csv here

            ###########################
            # process data

            logging.log("processing data")

            ta_sites_df,models_df,gene_obj_dict,filtered_ttn_data,gumbel_bernoulli_gene_calls = self.calc_ttnfitness(genome,G,self.inputs.gumbel_results_path)

            ###########################
            # write output
            # 
            # note: first header line is filetype, last header line is column headers

            self.write_ttnfitness_results(ta_sites_df,models_df,gene_obj_dict,filtered_ttn_data,gumbel_bernoulli_gene_calls,self.inputs.genes_output_path,self.inputs.sites_output_path) 


            if gui.is_active and self.inputs.genes_output_path!=None:
                logging.log(f"Adding File: {self.inputs.genes_output_path}")
                results_area.add(self.inputs.genes_output_path)
                logging.log(f"Adding File: {self.inputs.sites_output_path}")
                results_area.add(self.inputs.sites_output_path)
            logging.log("Finished TTNFitness")
            logging.log("Time: %0.1fs\n" % (time.time() - self.start_time))

    def calc_ttnfitness(self, genome, G, gumbel_results_file):
        """
        Returns:
            ta_sites_df, models_df, gene_obj_dict
        """
        import pandas
        import statsmodels.stats.multitest
        import statsmodels.api as sm
        
        self.gumbel_estimations = gumbel_results_file

        # Creating the dataset
        orf = []
        name = []
        coords = []
        ttn_vector_list = []
        # Get the nucleotides surrounding the TA sites
        genome2 = genome + genome
        all_counts = []
        combos = ["".join(p) for p in itertools.product(["A", "C", "T", "G"], repeat=4)]
        gene_obj_dict = {}
        upseq_list = []
        downseq_list = []
        for gene in G:
            gene_obj_dict[gene.orf] = gene
            all_counts.extend(
                numpy.mean(gene.reads, 0)
            )  # mean TA site counts across wig files
            for pos in gene.position:
                pos -= 1  # 1-based to 0-based indexing of nucleotides
                if pos - 4 < 0:
                    pos += len(genome)
                nucs = genome2[pos - 4 : pos + 6]
                if nucs[4:6] != "TA":
                    sys.stderr.write(
                        "warning: site %d is %s instead of TA" % (pos, nucs[4:6])
                    )
                # convert nucleotides to upstream and downstream TTN
                upseq = nucs[0] + nucs[1] + nucs[2] + nucs[3]
                upseq_list.append(upseq)
                # reverse complementing downstream
                downseq = ""
                for x in [nucs[9], nucs[8], nucs[7], nucs[6]]:
                    if str(x) == "A":
                        downseq += "T"
                    if str(x) == "C":
                        downseq += "G"
                    if str(x) == "G":
                        downseq += "C"
                    if str(x) == "T":
                        downseq += "A"
                downseq_list.append(downseq)
                ttn_vector = []
                for c in combos:
                    if upseq == c and downseq == c:
                        ttn_vector.append(int(2))  # up/dwn ttn are same, "bit"=2
                    elif upseq == c or downseq == c:
                        ttn_vector.append(int(1))  # set ttn bit=1
                    else:
                        ttn_vector.append(int(0))
                ttn_vector_list.append(pandas.Series(ttn_vector, index=combos))
                orf.append(gene.orf)
                name.append(gene.name)
                coords.append(pos)

        ta_sites_df = pandas.DataFrame(
            {
                "ORF": orf,
                "Gene Name": name,
                "Coordinates": coords,
                "Insertion Count": all_counts,
                "Upstream TTN": upseq_list,
                "Downstream TTN": downseq_list,
            }
        )
        ta_sites_df = pandas.concat(
            [ta_sites_df, pandas.DataFrame(ttn_vector_list)], axis=1
        )
        ta_sites_df = ta_sites_df.sort_values(by=["Coordinates"], ignore_index=True)
        # get initial states of the TA Sites
        # compute state labels (ES or NE)
        # for runs of >=R TA sites with cnt=0; label them as "ES", and the rest as "NE"
        # treat ends of genome as connected (circular)
        Nsites = len(ta_sites_df["Insertion Count"])
        states = ["NE"] * Nsites
        R = 6  # make this adaptive based on saturation?
        MinCount = 2
        i = 0
        while i < Nsites:
            j = i
            while j < Nsites and ta_sites_df["Insertion Count"].iloc[j] < MinCount:
                j += 1
            if j - i >= R:
                for k in range(i, j):
                    states[k] = "ES"
                i = j
            else:
                i += 1

        # getlocal averages -excludes self
        W = 5
        localmeans = []
        for i in range(Nsites):
            vals = []
            for j in range(-W, W + 1):
                if (
                    j != 0 and i + j >= 0 and i + j < Nsites
                ):  # this excludes the site itself
                    if states[i + j] != states[i]:
                        continue  # include only neighboring sites with same state when calculating localmean # diffs2.txt !!!
                    vals.append(float(ta_sites_df["Insertion Count"].iloc[i + j]))
            smoothed = -1 if len(vals) == 0 else numpy.mean(vals)
            localmeans.append(smoothed)

        # get LFCs
        lfc_values = []
        pseudocount = 10
        for i in range(len(ta_sites_df["Insertion Count"])):
            c, m = ta_sites_df["Insertion Count"].iloc[i], localmeans[i]
            lfc = math.log((c + pseudocount) / float(m + pseudocount), 2)
            lfc_values.append(lfc)

        ta_sites_df["State"] = states
        ta_sites_df["Local Average"] = localmeans
        ta_sites_df["Actual LFC"] = lfc_values

        ####################################################

        logging.log("Making Fitness Estimations")
        # Read in Gumbel estimations
        skip_count = 0
        gumbel_file = open(self.gumbel_estimations, "r")
        for line in gumbel_file.readlines():
            if line.startswith("#"):
                skip_count = skip_count + 1
            else:
                break
        gumbel_file.close()
        from .gumbel import Method as Gumbel
        gumbel_df = pandas.read_csv(
            self.gumbel_estimations,
            sep="\t",
            skiprows=skip_count,
            names=Gumbel.column_names,
            dtype=str,
        )

        saturation = len(ta_sites_df[ta_sites_df["Insertion Count"] > 0]) / len(ta_sites_df)
        phi = 1.0 - saturation
        significant_n = math.log10(0.05) / math.log10(phi)

        logging.log("\t + Filtering ES/ESB Genes")
        # function to extract gumbel calls to filter out ES and ESB
        gumbel_bernoulli_gene_calls = {}
        for _, g in informative_iterator.ProgressBar(ta_sites_df["ORF"].unique(), title="Filtering ES/ESB Genes"):
            if g == "igr":
                gene_call = numpy.nan
            else:
                gene_call = "U"
                sub_gumbel = gumbel_df[gumbel_df["ORF"] == g]
                if len(sub_gumbel) > 0:
                    gene_call = sub_gumbel["Essentiality Call"].iloc[0]
                # set to ES if greater than n and all 0s
                sub_data = ta_sites_df[ta_sites_df["ORF"] == g]
                if (
                    len(sub_data) > significant_n
                    and len(sub_data[sub_data["Insertion Count"] > 0]) == 0
                ):
                    gene_call = "EB"  # binomial filter
            gumbel_bernoulli_gene_calls[g] = gene_call
            percentage = _["percent"]
            if gui.is_active:
                parameter_panel.progress_update(f"Filtering ES/ESB Genes... {percentage:.2f}%", percentage)
        ess_genes = [
            key
            for key, value in gumbel_bernoulli_gene_calls.items()
            if (value == "E") or (value == "EB")
        ]
        
        logging.log("\t + Filtering Short Genes. Labeling as Uncertain")
        # function to call short genes (1 TA site) with no insertions as Uncertain
        uncertain_genes = []
        for _, g in informative_iterator.ProgressBar(ta_sites_df["ORF"].unique(), title="Filtering Short Genes"):
            sub_data = ta_sites_df[ta_sites_df["ORF"] == g]
            len_of_gene = len(sub_data)
            num_insertions = len(sub_data[sub_data["Insertion Count"] > 0])
            saturation = num_insertions / len_of_gene
            if saturation == 0 and len_of_gene <= 1:
                uncertain_genes.append(g)
            percentage = _["percent"]
            if gui.is_active:
                parameter_panel.progress_update(f"Filtering Short Genes. Labeling as Uncertain... {percentage:.2f}%", percentage)

        filtered_ttn_data = ta_sites_df[ta_sites_df["State"] != "ES"]
        filtered_ttn_data = filtered_ttn_data[filtered_ttn_data["Local Average"] != -1]
        filtered_ttn_data = filtered_ttn_data[
            ~filtered_ttn_data["ORF"].isin(ess_genes)
        ]  # filter out ess genes
        filtered_ttn_data = filtered_ttn_data[
            ~filtered_ttn_data["ORF"].isin(uncertain_genes)
        ]  # filter out uncertain genes
        filtered_ttn_data = filtered_ttn_data.reset_index(drop=True)


        ##########################################################################################
        # Linear Regression
        gene_one_hot_encoded = pandas.get_dummies(filtered_ttn_data["ORF"], prefix="").astype(int)
        #gene_one_hot_encoded.replace(to_replace = 0, value = -1/(len(gene_one_hot_encoded.columns)-1), inplace=True)

        columns_to_drop = [
            "Coordinates",
            "Insertion Count",
            "ORF",
            "Gene Name",
            "Local Average",
            "Upstream TTN",
            "Downstream TTN",
        ]
        assert all([ each in SitesFile.column_names for each in columns_to_drop]), f"Developer Error: TTN Fitness is dropping columns that were probably renamed: {[ each for each in columns_to_drop if each not in SitesFile.column_names]}"
        ttn_vectors = filtered_ttn_data.drop(
            [ "Actual LFC", "State",] + columns_to_drop,
            axis=1,
        )
        #ttn_vectors.replace(to_replace = 0, value = 1/(len(ttn_vectors.columns)-2), inplace=True)
   
        old_Y = numpy.log10(filtered_ttn_data["Insertion Count"] + 0.5)
        Y = old_Y - numpy.mean(old_Y) #centering Y values so we can disregard constant
        #Y = numpy.log10(filtered_ttn_data["Insertion Count"] + 0.5)

        
        percentage=50
        parameter_panel.progress_update(f"Fitting M1... {percentage:.2f}%", percentage)
        logging.log("\t + Fitting M1")
        if True: # NOTE: the block of code in this if statement is what takes up the bulk of the processing time (can't give good ETA/progress cause of this)
            X1, results1 = slow_computation(gene_one_hot_encoded, ttn_vectors, Y) # is a function so that it can be cached
            not_overlap =[]
            for i in gene_one_hot_encoded.columns:
                if i not in results1.params.index: not_overlap.append(i)
            not_overlap =[]
            for i in ttn_vectors.columns:
                if i not in results1.params.index: not_overlap.append(i)
            filtered_ttn_data["M1 Pred Log Count"] = results1.predict(X1) 
            filtered_ttn_data["M1 Pred Log Count"] = filtered_ttn_data["M1 Pred Log Count"] + numpy.mean(old_Y) #adding mean target value to account for centering
            filtered_ttn_data["M1 Predicted Count"] = numpy.power(
                10, (filtered_ttn_data["M1 Pred Log Count"] - 0.5)
            )

        logging.log("\t + Assessing Models")
        # create Models Summary df
        models_df = pandas.DataFrame(results1.params[1:-256], columns=["M1 Coef"])
        models_df["M1 P Value"] = results1.pvalues[1:-256]
        models_df["M1 Adj P Value"] = statsmodels.stats.multitest.fdrcorrection(
            results1.pvalues[1:-256], alpha=0.05
        )[1]

        # creating a mask for the adjusted pvals
        models_df.loc[
            (models_df["M1 Coef"] > 0) & (models_df["M1 Adj P Value"] < 0.05),
            "Gene Plus TTN States",
        ] = "GA"
        models_df.loc[
            (models_df["M1 Coef"] < 0) & (models_df["M1 Adj P Value"] < 0.05),
            "Gene Plus TTN States",
        ] = "GD"
        models_df.loc[
            (models_df["M1 Coef"] == 0) & (models_df["M1 Adj P Value"] < 0.05),
            "Gene Plus TTN States",
        ] = "NE"
        models_df.loc[(models_df["M1 Adj P Value"] > 0.05), "Gene Plus TTN States"] = "NE"

        return (ta_sites_df,models_df,gene_obj_dict,filtered_ttn_data,gumbel_bernoulli_gene_calls)

    def write_ttnfitness_results(self, ta_sites_df, models_df, gene_obj_dict, filtered_ttn_data, gumbel_bernoulli_gene_calls, genes_output_path, sites_output_path):
        import pandas
        genes_out_rows, sites_out_rows = [],[]
        logging.log("Writing To Output Files")
        # Write Models Information to CSV
        # Columns: ORF ID, ORF Name, ORF Description,M0 Coef, M0 Adj P Value

        gene_dict = {}  # dictionary to map information per gene
        ta_sites_df["M1 Predicted Count"] = [None] * len(ta_sites_df)
        for progress, g in informative_iterator.ProgressBar(ta_sites_df["ORF"].unique(), title="Writing To Output"):
            # Update progress
            percentage = progress["percent"]
            if gui.is_active:
                parameter_panel.progress_update(f"Writing Output Files... {percentage:.2f}%", percentage)
            # ORF Name
            orf_name = gene_obj_dict[g].name
            # ORF Description
            orf_description = gene_obj_dict[g].desc
            # Total TA sites
            num_t_asites = len(gene_obj_dict[g].reads[0])  # TRI check this!
            # Sites > 0
            above0_ta_sites = len([r for r in gene_obj_dict[g].reads[0] if r > 0])
            # Insertion Count
            actual_counts = ta_sites_df[ta_sites_df["ORF"] == g]["Insertion Count"]
            mean_actual_counts = numpy.mean(actual_counts)
            local_saturation = above0_ta_sites / num_t_asites
            # Predicted Count
            if g in filtered_ttn_data["ORF"].values:
                actual_df = filtered_ttn_data[filtered_ttn_data["ORF"] == g][
                    "Insertion Count"
                ]
                coords_orf = filtered_ttn_data[filtered_ttn_data["ORF"] == g][
                    "Coordinates"
                ].values.tolist()
                for c in coords_orf:
                    ta_sites_df.loc[
                        (ta_sites_df["Coordinates"] == c), "M1 Predicted Count"
                    ] = filtered_ttn_data[filtered_ttn_data["Coordinates"] == c][
                        "M1 Predicted Count"
                    ].iloc[
                        0
                    ]
            # M1 info
            if "_" + g in models_df.index:
                m1_coef = models_df.loc["_" + g, "M1 Coef"]
                m1_adj_p_value = models_df.loc["_" + g, "M1 Adj P Value"]
                modified_m1 = math.exp(
                    m1_coef - statistics.median(models_df["M1 Coef"].values.tolist())
                )
            else:
                m1_coef = ""
                m1_adj_p_value = ""
                modified_m1 = ""

            # States
            gumbel_bernoulli_call = gumbel_bernoulli_gene_calls[g]
            if gumbel_bernoulli_call == "E":
                gene_ttn_call = "ES"
            elif gumbel_bernoulli_call == "EB":
                gene_ttn_call = "ESB"
            else:
                if "_" + g in models_df.index:
                    gene_ttn_call = models_df.loc["_" + g, "Gene Plus TTN States"]
                else:
                    gene_ttn_call = "U"  # these genes are in the uncertain genes list
            ta_sites_df.loc[
                (ta_sites_df["ORF"] == g), "TTN Fitness Assessment"
            ] = gene_ttn_call
            gene_dict[g] = [
                g,
                orf_name,
                orf_description,
                num_t_asites,
                above0_ta_sites,
                local_saturation,
                m1_coef,
                m1_adj_p_value,
                mean_actual_counts,
                modified_m1,
                gene_ttn_call,
            ]
        

        saturation = len(ta_sites_df[ta_sites_df["Insertion Count"] > 0]) / len(ta_sites_df) 
        

        
        # 
        # Write Genes data
        # 
        output_df = pandas.DataFrame.from_dict(gene_dict, orient="index")
        output_df.columns = GenesFile.column_names
        assesment_cnt = output_df["TTN Fitness Assessment"].value_counts()
        genes_out_rows = output_df.values.tolist()
        
        #in case GA, GD etc. do not exist
        summary_info={}
        for ess in ["ES", "ESB","GD","GA","NE","U"]:
            if ess not in assesment_cnt: summary_info[ess]="0"
            else: summary_info[ess] = str(assesment_cnt[ess])

        logging.log("Writing File: %s" % (self.inputs.genes_output_path))
        transit_tools.write_result(
            path=Method.inputs.genes_output_path,
            file_kind=Method.identifier+"Genes",
            rows=genes_out_rows,
            column_names=output_df.columns,
            extra_info=dict(
                condition_selected = Method.inputs.condition,
                calculation_time=f"{(time.time() - self.start_time):0.1f}seconds",
                saturation= saturation,
                analysis_type=Method.identifier,
                files=dict(
                    combined_wig = Method.inputs.combined_wig,
                    wig_files = Method.inputs.wig_files,
                    metadata = Method.inputs.metadata,
                    annotation_path=Method.inputs.annotation_path,
                    gumbel_results_file = Method.inputs.gumbel_results_path,
                ),
                parameters=dict(
                    normalization = self.inputs.normalization,
                ),
                summary_info=summary_info,

                naming_reference = dict(
                    ES = "Essential based on Gumbel",
                    ESB = "Essential based on Binomial",
                    GD = "Growth Defect",
                    NE = "Non-essential",
                    GA = "Growth Advantage",
                    U = "Uncertain",
                ),
            ),
        )
        
        # 
        # write Sites data
        # 
        ta_sites_df = ta_sites_df[SitesFile.column_names]
        sites_out_rows = ta_sites_df.values.tolist()
        logging.log("Writing File: %s" % (self.inputs.sites_output_path))
        transit_tools.write_result(
            path=self.inputs.sites_output_path,
            file_kind=Method.identifier+"Sites",
            rows=sites_out_rows,
            column_names=ta_sites_df.columns,
            extra_info=dict(
                condition_selected = Method.inputs.condition,
                calculation_time=f"{(time.time() - self.start_time):0.1f}seconds",
                analysis_type=Method.identifier,
                saturation = saturation,
                files=dict(
                    combined_wig = self.inputs.combined_wig,
                    wig_files = self.inputs.wig_files,
                    metadata = self.inputs.metadata,
                    annotation_path=self.inputs.annotation_path,
                    gumbel_results_file = self.inputs.gumbel_results_path,
                ),
                parameters=dict(
                    normalization = self.inputs.normalization,
                ),
            ),
        )

        logging.log(f"Finished running {Method.identifier}")  # Printing empty line to flush stdout


            

@transit_tools.ResultsFile
class GenesFile:
    column_names = [
        "ORF",
        "Gene Name",
        "Description",
        "Total TA Site Count",
        "Count Of Sites With Insertions",
        "Gene Saturation",
        "Gene Plus TTN M1 Coef",
        "Gene Plus TTN M1 Adj P Value",
        "Mean Insertion Count",
        "Fitness Ratio",
        "TTN Fitness Assessment",
    ]
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier+"Genes")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=transit_tools.basename(self.path),
            type=Method.identifier+"Genes",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier,
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                ).Show(),

                "Display Volcano Plot": lambda *args: self.graph_volcano_plot(),
            })
        )
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        
        summary = self.extra_data.get("summary_info", {})
        summary_str = [str(summary[key])+" "+str(key) for key in ["ES", "ESB", "GD", "GA", "NE", "U"]] 
        self.values_for_result_table.update({"summary": "; ".join(summary_str) })
        
        parameters = self.extra_data.get("parameters",{})
        parameters_str = [str(key)+" : "+str(parameters[key]) for key in ["normalization"]]
        self.values_for_result_table.update({"parameters": "; ".join(parameters_str) })

    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()


    def graph_volcano_plot(self):
        import pandas
        with gui_tools.nice_error_log:
            Method.inputs.volcano_output_path = gui_tools.ask_for_output_file_path(
                default_file_name=f"volcano.png",
                output_extensions='PNG file (*.png)|*.png;|\nAll files (*.*)|*.*',
            )

            try: import matplotlib.pyplot as plt
            except:
                print("Error: cannot do plots, no matplotlib")

            ttnfitness_genes_summary = pandas.read_csv(self.path, sep= "\t",comment='#')
            ttnfitness_genes_summary.columns = GenesFile.column_names

            color_dict = {
                "ES":"r",
                "ESB" : "b",
                "NE" : "g",
                "GA" : "m",
                "GD" : "c",
                "U": "y"
            }
            plt.figure()
            for call in ["ES","ESB", "GD","GA", "NE","U"]:
                sub_summary = ttnfitness_genes_summary[ttnfitness_genes_summary["TTN Fitness Assessment"]==call]
                coef_vals = sub_summary["Gene Plus TTN M1 Coef"]
                q_vals = sub_summary["Gene Plus TTN M1 Adj P Value"]
                log10_q_vals = []
                for each_q_val in q_vals:
                    try:
                        log10_q_value = -math.log(float(each_q_val), 10)
                    except ValueError as e:
                        log10_q_value = None
                    
                    log10_q_vals.append(log10_q_value)
                threshold = 0.05

                plt.scatter(coef_vals, log10_q_vals, c = color_dict[call], marker=".", label=call, alpha=0.5)
            plt.axhline( -math.log(threshold, 10), color="k", linestyle="dashed", linewidth=2)
            plt.axvline(0, color="k", linestyle="dashed", linewidth=2)
            plt.legend()
            plt.xlabel("Gene Plus TTN M1 Coef")
            plt.ylabel("-Log Adj P Value (base 10)")
            #plt.suptitle("Volcano plot")
            plt.title("Adjusted threshold (horizonal line): P-value=%1.8f\nVertical line set at Coef=0" % threshold)
            #plt.show()
            plt.savefig(Method.inputs.volcano_output_path, bbox_inches='tight')

            if gui.is_active:
                logging.log(f"Adding File: {Method.inputs.volcano_output_path}")
                results_area.add(Method.inputs.volcano_output_path)
            

@transit_tools.ResultsFile
class SitesFile:
    column_names = [
        "Coordinates",
        "ORF",
        "Gene Name",
        "Upstream TTN",
        "Downstream TTN",
        "TTN Fitness Assessment",
        "Insertion Count",
        "Local Average",
        "M1 Predicted Count",
    ]
    
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier+"Sites")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=transit_tools.basename(self.path),
            type=Method.identifier+"Sites",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier,
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                ).Show(),
            })
        )

        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
