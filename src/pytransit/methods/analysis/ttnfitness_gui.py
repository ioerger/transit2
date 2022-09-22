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

import pandas
import statsmodels.stats.multitest
import statsmodels.api as sm

from pytransit.basics.lazy_dict import LazyDict

from pytransit.universal_data import universal
from pytransit.components.parameter_panel import panel as parameter_panel
from pytransit.components.parameter_panel import progress_update
from pytransit.components.panel_helpers import *
from pytransit.components.spreadsheet import SpreadSheet
import pytransit.tools.gui_tools as gui_tools
import pytransit.tools.transit_tools as transit_tools
import pytransit.tools.tnseq_tools as tnseq_tools
import pytransit.tools.norm_tools as norm_tools
import pytransit.tools.stat_tools as stat_tools
import pytransit.basics.csv as csv
import pytransit.components.results_area as results_area

command_name = sys.argv[0]

class Analysis:
    identifier  = "TTNFitness"
    short_name  = "ttnfitness_gui"
    long_name   = "ttnfitness_gui"
    short_desc  = "Analyze fitness effect of (non-essential) genes"
    long_desc   = """Analyze fitness effect of (non-essential) genes using a predictive model that corrects for the bias in Himar1 insertion probabilities based on nucleotides around each TA site"""
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
    

    usage_string = f"""usage: python3 %s ttnfitness_gui <comma-separated .wig files> <annotation .prot_table> <genome .fna> <gumbel results file> <genes output file> <sites output file>""" % sys.argv[0] # the old way, with multiple wigs as input
    #usage_string = f"""usage: python3 %s ttnfitness <combined_wig> <sample_metadata> <condition> <gumbel_output_file> <genome .fna> <annotation .prot_table> <gumbel output file> <genes_output_file> <sites_output_file>""" # add '-c' to indicate combined_wig?
    
    wxobj = None
    panel = None
    
    def __init__(self, *args, **kwargs):
        Analysis.instance = self
        self.full_name        = f"[{self.short_name}]  -  {self.short_desc}"
        self.transposons_text = transit_tools.get_transposons_text(self.transposons)
        self.filetypes        = [GenesFile, SitesFile]
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

    def create_input_field(self, panel, sizer, label, value,tooltip=None):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text=label,
            default_value=value,
            tooltip_text=tooltip,
        )
        return lambda *args: get_text()

    def define_panel(self, _):
        from pytransit.tools.transit_tools import wx
        self.panel = make_panel()
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.value_getters = LazyDict()

        self.value_getters.condition = create_condition_choice(self.panel,main_sizer,"Condition to analyze:")
        self.value_getters.gumbel_results_path = create_file_input(self.panel,main_sizer, \
          button_label="Gumbel results file",default_file_name="glycerol_gumbel.out",allowed_extensions="All files (*.*)|*.*", \
          popup_title="Choose Gumbel results file", \
          tooltip_text="Must run Gumbel first to determine which genes are essential. Note: TTN-fitness estimates fitness of NON-essential genes.")
        self.value_getters.genome_path = create_file_input(self.panel,main_sizer, \
          popup_title="Choose genome sequence file", \
          button_label="Load genome sequence file",default_file_name="H37Rv.fna",allowed_extensions="Fasta files (*.fa;*.fna;*.fasta))|*.fa;*.fna;*.fasta", \
          tooltip_text="Genome sequence file (.fna) must match annotation file (.prot_table)")
        self.value_getters.output_basename = self.create_input_field(self.panel,main_sizer, \
          label="Basename for output files",value="ttnfitness.test",tooltip="If X is basename, then X_genes.dat and X_sites.dat will be generated as output files.")
        self.value_getters.normalization = create_normalization_input(self.panel, main_sizer) # TTR 

        create_run_button(self.panel, main_sizer)

        parameter_panel.set_panel(self.panel)
        self.panel.SetSizer(main_sizer)
        self.panel.Layout()
        main_sizer.Fit(self.panel)

    @classmethod
    def from_gui(cls, frame):
        with gui_tools.nice_error_log:
            combined_wig = universal.session_data.combined_wigs[0]
            Analysis.inputs.combined_wig = combined_wig.main_path
            # assume all samples are in the same metadata file
            Analysis.inputs.metadata_path = universal.session_data.combined_wigs[0].metadata_path 

            Analysis.inputs.annotation_path = universal.session_data.annotation_path

            for each_key, each_getter in Analysis.instance.value_getters.items():
                try:
                    Analysis.inputs[each_key] = each_getter()
                except Exception as error:
                    raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
            ###logging.log("included_conditions", Analysis.inputs.included_conditions)
            Analysis.inputs.genes_output_path = "%s.genes.dat" % (Analysis.inputs.output_basename)
            Analysis.inputs.sites_output_path = "%s.sites.dat" % (Analysis.inputs.output_basename)


            #if not Analysis.inputs.output_path: return None ### why?
            return Analysis.instance

    @classmethod
    def from_args(cls, args, kwargs): # clean_args() was already called in pytransit/__main__.py

      if len(args)!=6: print(cls.usage_string); sys.exit(0) # use transit_error()?

      Analysis.inputs.update(dict(
        combined_wig = None,
        metadata = None,
        wig_files = args[0].split(','),
        annotation_path = args[1],
        genome_path = args[2],
        gumbel_results_path = args[3],
        genes_output_path = args[4],
        sites_output_path = args[5],
        #normalization = "TTR",
      ))
        
      return Analysis.instance
        
    def Run(self):
        with gui_tools.nice_error_log:
            logging.log("Starting tnseq_stats analysis")
            self.start_time = time.time()

            #######################
            # get data

            if self.inputs.combined_wig!=None:  # assume metadata and condition are defined too
              logging.log("Getting Data from %s" % self.inputs.combined_wig)
              position, data, filenames_in_comb_wig = tnseq_tools.read_combined_wig(self.inputs.combined_wig)

              metadata = tnseq_tools.CombinedWigMetadata(self.inputs.metadata_path)
              indexes = {}
              for i,row in enumerate(metadata.rows): 
                cond = row["Condition"] 
                if cond not in indexes: indexes[cond] = []
                indexes[cond].append(i)
              cond = Analysis.inputs.condition
              ids = [metadata.rows[i]["Id"] for i in indexes[cond]]
              logging.log("selected samples for ttnfitness (cond=%s): %s" % (cond,','.join(ids)))
              data = data[indexes[cond]] # project array down to samples selected by condition

              # now, select the columns in data corresponding to samples that are replicates of desired condition...

            elif self.inputs.wig_files!=None:
              logging.log("Getting Data")
              (data, position) = transit_tools.get_validated_data( self.inputs.wig_files, wxobj=self.wxobj )

            else: print("error: must provide either combined_wig or list of wig files"); sys.exit(0) ##### use transit.error()?
                
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
                #minread=1,  ### add these options?
                #reps=self.replicates,
                #ignore_codon=self.ignore_codon,
                #n_terminus=self.n_terminus, 
                #c_terminus=self.c_terminus,
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

            TA_sites_df,Models_df,gene_obj_dict,filtered_ttn_data,gumbel_bernoulli_gene_calls = self.calc_ttnfitness(genome,G,self.inputs.gumbel_results_path)

            ###########################
            # write output
            # 
            # note: first header line is filetype, last header line is column headers

            self.write_ttnfitness_results(TA_sites_df,Models_df,gene_obj_dict,filtered_ttn_data,gumbel_bernoulli_gene_calls,self.inputs.genes_output_path,self.inputs.sites_output_path) 


            if universal.interface=="gui" and self.inputs.genes_output_path!=None:
              logging.log(f"Adding File: {self.inputs.genes_output_path}")
              results_area.add(self.inputs.genes_output_path)
              logging.log(f"Adding File: {self.inputs.sites_output_path}")
              results_area.add(self.inputs.sites_output_path)
            logging.log("Finished TnseqStats")
            logging.log("Time: %0.1fs\n" % (time.time() - self.start_time))

    # returns: TA_sites_df , Models_df , gene_obj_dict

    def calc_ttnfitness(self,genome,G,gumbel_results_file):
        self.gumbelestimations = gumbel_results_file

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
            nTAs = len(gene.reads[0])
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

        TA_sites_df = pandas.DataFrame(
            {
                "Orf": orf,
                "Name": name,
                "Coord": coords,
                "Insertion Count": all_counts,
                "Upstream TTN": upseq_list,
                "Downstream TTN": downseq_list,
            }
        )
        TA_sites_df = pandas.concat(
            [TA_sites_df, pandas.DataFrame(ttn_vector_list)], axis=1
        )
        TA_sites_df = TA_sites_df.sort_values(by=["Coord"], ignore_index=True)
        # get initial states of the TA Sites
        # compute state labels (ES or NE)
        # for runs of >=R TA sites with cnt=0; label them as "ES", and the rest as "NE"
        # treat ends of genome as connected (circular)
        Nsites = len(TA_sites_df["Insertion Count"])
        states = ["NE"] * Nsites
        R = 6  # make this adaptive based on saturation?
        MinCount = 2
        i = 0
        while i < Nsites:
            j = i
            while j < Nsites and TA_sites_df["Insertion Count"].iloc[j] < MinCount:
                j += 1
            if j - i >= R:
                for k in range(i, j):
                    states[k] = "ES"
                i = j
            else:
                i += 1

        # getlocal averages --excludes self
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
                    vals.append(float(TA_sites_df["Insertion Count"].iloc[i + j]))
            smoothed = -1 if len(vals) == 0 else numpy.mean(vals)
            localmeans.append(smoothed)

        # get LFCs
        LFC_values = []
        pseudocount = 10
        for i in range(len(TA_sites_df["Insertion Count"])):
            c, m = TA_sites_df["Insertion Count"].iloc[i], localmeans[i]
            lfc = math.log((c + pseudocount) / float(m + pseudocount), 2)
            LFC_values.append(lfc)

        TA_sites_df["State"] = states
        TA_sites_df["Local Average"] = localmeans
        TA_sites_df["Actual LFC"] = LFC_values

        ####################################################

        logging.log("Making Fitness Estimations")
        # Read in Gumbel estimations
        skip_count = 0
        gumbel_file = open(self.gumbelestimations, "r")
        for line in gumbel_file.readlines():
            if line.startswith("#"):
                skip_count = skip_count + 1
            else:
                break
        gumbel_file.close()
        gumbel_df = pandas.read_csv(
            self.gumbelestimations,
            sep="\t",
            skiprows=skip_count,
            names=["Orf", "Name", "Desc", "k", "n", "r", "s", "zbar", "Call"],
            dtype=str,
        )

        saturation = len(TA_sites_df[TA_sites_df["Insertion Count"] > 0]) / len(TA_sites_df)
        phi = 1.0 - saturation
        significant_n = math.log10(0.05) / math.log10(phi)

        logging.log("\t + Filtering ES/ESB Genes")
        # function to extract gumbel calls to filter out ES and ESB
        gumbel_bernoulli_gene_calls = {}
        for g in TA_sites_df["Orf"].unique():
            if g == "igr":
                gene_call = numpy.nan
            else:
                gene_call = "U"
                sub_gumbel = gumbel_df[gumbel_df["Orf"] == g]
                if len(sub_gumbel) > 0:
                    gene_call = sub_gumbel["Call"].iloc[0]
                # set to ES if greater than n and all 0s
                sub_data = TA_sites_df[TA_sites_df["Orf"] == g]
                if (
                    len(sub_data) > significant_n
                    and len(sub_data[sub_data["Insertion Count"] > 0]) == 0
                ):
                    gene_call = "EB"  # binomial filter
            gumbel_bernoulli_gene_calls[g] = gene_call
        ess_genes = [
            key
            for key, value in gumbel_bernoulli_gene_calls.items()
            if (value == "E") or (value == "EB")
        ]

        logging.log("\t + Filtering Short Genes. Labeling as Uncertain")
        # function to call short genes (1 TA site) with no insertions as Uncertain
        uncertain_genes = []
        for g in TA_sites_df["Orf"].unique():
            sub_data = TA_sites_df[TA_sites_df["Orf"] == g]
            len_of_gene = len(sub_data)
            num_insertions = len(sub_data[sub_data["Insertion Count"] > 0])
            saturation = num_insertions / len_of_gene
            if saturation == 0 and len_of_gene <= 1:
                uncertain_genes.append(g)

        filtered_ttn_data = TA_sites_df[TA_sites_df["State"] != "ES"]
        filtered_ttn_data = filtered_ttn_data[filtered_ttn_data["Local Average"] != -1]
        filtered_ttn_data = filtered_ttn_data[
            ~filtered_ttn_data["Orf"].isin(ess_genes)
        ]  # filter out ess genes
        filtered_ttn_data = filtered_ttn_data[
            ~filtered_ttn_data["Orf"].isin(uncertain_genes)
        ]  # filter out uncertain genes
        filtered_ttn_data = filtered_ttn_data.reset_index(drop=True)


        ##########################################################################################
        # Linear Regression
        gene_one_hot_encoded = pandas.get_dummies(filtered_ttn_data["Orf"], prefix="")
        ttn_vectors = filtered_ttn_data.drop(
            [
                "Coord",
                "Insertion Count",
                "Orf",
                "Name",
                "Local Average",
                "Actual LFC",
                "State",
                "Upstream TTN",
                "Downstream TTN",
            ],
            axis=1,
        )
   
        old_Y = numpy.log10(filtered_ttn_data["Insertion Count"] + 0.5)
        Y = old_Y - numpy.mean(old_Y) #centering Y values so we can disregard constant



        logging.log("\t + Fitting M1")
        X1 = pandas.concat([gene_one_hot_encoded, ttn_vectors], axis=1)
        #X1 = sm.add_constant(X1)
        results1 = sm.OLS(Y, X1).fit()
        filtered_ttn_data["M1 Pred log Count"] = results1.predict(X1) 
        filtered_ttn_data["M1 Pred log Count"] = filtered_ttn_data["M1 Pred log Count"] + numpy.mean(old_Y) #adding mean target value to account for centering
        filtered_ttn_data["M1 Predicted Count"] = numpy.power(
            10, (filtered_ttn_data["M1 Pred log Count"] - 0.5)
        )

        logging.log("\t + Assessing Models")
        # create Models Summary df
        Models_df = pandas.DataFrame(results1.params[1:-256], columns=["M1 Coef"])
        Models_df["M1 Pval"] = results1.pvalues[1:-256]
        Models_df["M1 Adjusted Pval"] = statsmodels.stats.multitest.fdrcorrection(
            results1.pvalues[1:-256], alpha=0.05
        )[1]

        # creating a mask for the adjusted pvals
        Models_df.loc[
            (Models_df["M1 Coef"] > 0) & (Models_df["M1 Adjusted Pval"] < 0.05),
            "Gene+TTN States",
        ] = "GA"
        Models_df.loc[
            (Models_df["M1 Coef"] < 0) & (Models_df["M1 Adjusted Pval"] < 0.05),
            "Gene+TTN States",
        ] = "GD"
        Models_df.loc[
            (Models_df["M1 Coef"] == 0) & (Models_df["M1 Adjusted Pval"] < 0.05),
            "Gene+TTN States",
        ] = "NE"
        Models_df.loc[(Models_df["M1 Adjusted Pval"] > 0.05), "Gene+TTN States"] = "NE"

        return (TA_sites_df,Models_df,gene_obj_dict,filtered_ttn_data,gumbel_bernoulli_gene_calls)

    def write_ttnfitness_results(self,TA_sites_df,Models_df,gene_obj_dict,filtered_ttn_data,gumbel_bernoulli_gene_calls,genes_output_path,sites_output_path):
        genes_out_rows, sites_out_rows = [],[]
        logging.log("Writing To Output Files")
        # Write Models Information to CSV
        # Columns: ORF ID, ORF Name, ORF Description,M0 Coef, M0 Adj Pval

        gene_dict = {}  # dictionary to map information per gene
        TA_sites_df["M1 Predicted Count"] = [None] * len(TA_sites_df)
        # TA_sites_df["mod ttn Predicted Count"] = [None]*len(TA_sites_df)
        for g in TA_sites_df["Orf"].unique():
            # ORF Name
            orfName = gene_obj_dict[g].name
            # ORF Description
            orfDescription = gene_obj_dict[g].desc
            # Total TA sites
            numTAsites = len(gene_obj_dict[g].reads[0])  # TRI check this!
            # Sites > 0
            above0TAsites = len([r for r in gene_obj_dict[g].reads[0] if r > 0])
            # Insertion Count
            actual_counts = TA_sites_df[TA_sites_df["Orf"] == g]["Insertion Count"]
            mean_actual_counts = numpy.mean(actual_counts)
            local_saturation = above0TAsites / numTAsites
            # Predicted Count
            if g in filtered_ttn_data["Orf"].values:
                actual_df = filtered_ttn_data[filtered_ttn_data["Orf"] == g][
                    "Insertion Count"
                ]
                coords_orf = filtered_ttn_data[filtered_ttn_data["Orf"] == g][
                    "Coord"
                ].values.tolist()
                for c in coords_orf:
                    TA_sites_df.loc[
                        (TA_sites_df["Coord"] == c), "M1 Predicted Count"
                    ] = filtered_ttn_data[filtered_ttn_data["Coord"] == c][
                        "M1 Predicted Count"
                    ].iloc[
                        0
                    ]
            # M1 info
            if "_" + g in Models_df.index:
                M1_coef = Models_df.loc["_" + g, "M1 Coef"]
                M1_adj_pval = Models_df.loc["_" + g, "M1 Adjusted Pval"]
                modified_M1 = math.exp(
                    M1_coef - statistics.median(Models_df["M1 Coef"].values.tolist())
                )
            else:
                M1_coef = None
                M1_adj_pval = None
                modified_M1 = None

            # States
            gumbel_bernoulli_call = gumbel_bernoulli_gene_calls[g]
            if gumbel_bernoulli_call == "E":
                gene_ttn_call = "ES"
            elif gumbel_bernoulli_call == "EB":
                gene_ttn_call = "ESB"
            else:
                if "_" + g in Models_df.index:
                    gene_ttn_call = Models_df.loc["_" + g, "Gene+TTN States"]
                else:
                    gene_ttn_call = "U"  # these genes are in the uncertain genes list
            TA_sites_df.loc[
                (TA_sites_df["Orf"] == g), "TTN-Fitness Assessment"
            ] = gene_ttn_call
            gene_dict[g] = [
                g,
                orfName,
                orfDescription,
                numTAsites,
                above0TAsites,
                local_saturation,
                M1_coef,
                M1_adj_pval,
                mean_actual_counts,
                modified_M1,
                gene_ttn_call,
            ]
        output_df = pandas.DataFrame.from_dict(gene_dict, orient="index")
        output_df.columns = [
            "ORF ID",
            "Name",
            "Description",
            "Total # TA Sites",
            "#Sites with insertions",
            "Gene Saturation",
            "Gene+TTN (M1) Coef",
            "Gene+TTN (M1) Adj Pval",
            "Mean Insertion Count",
            "Fitness Ratio",
            "TTN-Fitness Assessment",
        ]
        assesment_cnt = output_df["TTN-Fitness Assessment"].value_counts()

        saturation = len(TA_sites_df[TA_sites_df["Insertion Count"] > 0]) / len(TA_sites_df) 
        TA_sites_df = TA_sites_df[
            [
                "Coord",
                "Orf",
                "Name",
                "Upstream TTN",
                "Downstream TTN",
                "TTN-Fitness Assessment",
                "Insertion Count",
                "Local Average",
                "M1 Predicted Count",
            ]
        ]

        genes_out_rows = output_df.values.tolist()

        logging.log("Writing File: %s" % (self.inputs.genes_output_path))
        transit_tools.write_result(
            path=self.inputs.genes_output_path,
            file_kind=Analysis.identifier+"Genes",
            rows=genes_out_rows,
            column_names=output_df.columns,
            extra_info=dict(
                parameters=dict(
                    combined_wig = self.inputs.combined_wig,
                    wig_files = self.inputs.wig_files,
                    metadata = self.inputs.metadata,
                    annotation_path=self.inputs.annotation_path,
                    gumbel_results_file = self.inputs.gumbel_results_path,
                    normalization = self.inputs.normalization,
                ),
                time=(time.time() - self.start_time),
                saturation = saturation,

                ES = str(assesment_cnt["ES"]) + " #essential based on Gumbel",
                ESB = str(assesment_cnt["ESB"]) + " #essential based on Binomial",
                GD = str(assesment_cnt["GD"]) +" #Growth Defect",
                GA = str(assesment_cnt["GA"]) +" #Growth Advantage",
                NE = str(assesment_cnt["NE"]) + " #non-essential",
                U = str(assesment_cnt["U"]) + " #uncertain",
                
                
            ),
        )
        # write sites data

        sites_out_rows = TA_sites_df.values.tolist()

        logging.log("Writing File: %s" % (self.inputs.sites_output_path))
        transit_tools.write_result(
            path=self.inputs.sites_output_path,
            file_kind=Analysis.identifier+"Sites",
            rows=sites_out_rows,
            column_names=TA_sites_df.columns,
            extra_info=dict(
                parameters=dict(
                    combined_wig = self.inputs.combined_wig,
                    wig_files = self.inputs.wig_files,
                    metadata = self.inputs.metadata,
                    annotation_path=self.inputs.annotation_path,
                    gumbel_results_file = self.inputs.gumbel_results_path,
                    normalization = self.inputs.normalization,
                ),
                time=(time.time() - self.start_time),
                saturation = saturation,

                ES = str(assesment_cnt["ES"]) + " #essential based on Gumbel",
                ESB = str(assesment_cnt["ESB"]) + " #essential based on Binomial",
                GD = str(assesment_cnt["GD"]) +" #Growth Defect",
                GA = str(assesment_cnt["GA"]) +" #Growth Advantage",
                NE = str(assesment_cnt["NE"]) + " #non-essential",
                U = str(assesment_cnt["U"]) + " #uncertain",        
            ),
        )

        logging.log("")  # Printing empty line to flush stdout
        # logging.log("Adding File: %s" % (self.output.name))
        # results_area.add(self.output.name)
        #self.finish()
        logging.log("Finished TTNFitness Method")


            

@transit_tools.ResultsFile
class GenesFile(Analysis):
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Analysis.identifier+"Genes")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=transit_tools.basename(self.path),
            type=Analysis.identifier+"Genes",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title="TTNFitness Summary",heading="",column_names=self.column_names,rows=self.rows).Show(),
                "Display Volcano Plot": lambda *args: self.graph_volcano_plot(),
            })
        )
        
        # 
        # get column names
        # 
        comments, headers, rows = csv.read(self.path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
        if len(comments) == 0:
            raise Exception(f'''No comments in file, and I expected the last comment to be the column names, while to load tnseq_stats file "{self.path}"''')
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


    def graph_volcano_plot(self):
        with gui_tools.nice_error_log:
            try: import matplotlib.pyplot as plt
            except:
                print("Error: cannot do plots, no matplotlib")

            ttnfitness_genes_summary = pandas.read_csv(self.path, sep= "\t",comment='#')
            ttnfitness_genes_summary.columns = [
                "ORF ID",
                "Name",
                "Description",
                "Total # TA Sites",
                "#Sites with insertions",
                "Gene Saturation",
                "Gene+TTN (M1) Coef",
                "Gene+TTN (M1) Adj Pval",
                "Mean Insertion Count",
                "Fitness Ratio",
                "TTN-Fitness Assessment",
            ]

            color_dict = {
                "ES":"r",
                "ESB" : "b",
                "NE" : "g",
                "GA" : "m",
                "GD" : "c",
                "U": "y"
            }
            plt.figure()
            for call in set(ttnfitness_genes_summary["TTN-Fitness Assessment"]):
                sub_summary = ttnfitness_genes_summary[ttnfitness_genes_summary["TTN-Fitness Assessment"]==call]
                coef_vals = sub_summary["Gene+TTN (M1) Coef"]
                q_vals = sub_summary["Gene+TTN (M1) Adj Pval"]
                log10_q_vals = []
                for each_q_val in q_vals:
                    try:
                        log10_q_value = -math.log(float(each_q_val), 10)
                    except ValueError as e:
                        log10_q_value = None
                    
                    log10_q_vals.append(log10_q_value)
                threshold = 0.05

                plt.scatter(coef_vals, log10_q_vals, c = color_dict[call], marker=".", label=call)
            plt.axhline( -math.log(threshold, 10), color="k", linestyle="dashed", linewidth=2)
            plt.axvline(0, color="k", linestyle="dashed", linewidth=2)
            plt.legend()
            plt.xlabel("Gene+TTN (M1) Coef")
            plt.ylabel("-Log adjusted p-value (base 10)")
            plt.suptitle("Resampling - Volcano plot")
            plt.title("Adjusted threshold (horizonal line): P-value=%1.8f\nVertical line set at Coef=0" % threshold)
            plt.show()
            

@transit_tools.ResultsFile
class SitesFile(Analysis):
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Analysis.identifier+"Sites")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=transit_tools.basename(self.path),
            type=Analysis.identifier+"Sites",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title="TTNFitness Summary",heading="",column_names=self.column_names,rows=self.rows).Show(),
            })
        )
        
        # 
        # get column names
        # 
        comments, headers, rows = csv.read(self.path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
        if len(comments) == 0:
            raise Exception(f'''No comments in file, and I expected the last comment to be the column names, while to load tnseq_stats file "{self.path}"''')
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
