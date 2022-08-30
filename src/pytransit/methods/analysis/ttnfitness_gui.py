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
    identifier  = "#ttnfitness_gui"
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
    
    valid_cli_flags = [
    # -n for normalization?
    ]

    usage_string = f"""usage: python3 %s ttnfitness <comma-separated .wig files> <annotation .prot_table> <genome .fna> <gumbel output file> <output1 file> <output2 file>""" # the old way, with multiple wigs as input
    #usage_string = f"""usage: python3 %s ttnfitness <combined_wig> <sample_metadata> <condition> <gumbel_output_file> <genome .fna> <annotation .prot_table> <gumbel output file> <genes_output_file> <sites_output_file>"""
    
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

        self.value_getters.gumbel_results_path = self.create_input_field(self.panel,main_sizer,label="Gumbel results file:",value="",tooltip="Must run Gumbel first to determine which genes are essential. Note: TTN-fitness estimates fitness of NON-essential genes.")
        self.value_getters.genome_path = self.create_input_field(self.panel,main_sizer,label="Genome sequence file:",value="",tooltip="For example, a .fasta or .fna file.")
        self.value_getters.output_basename = self.create_input_field(self.panel,main_sizer,label="Basename for output files",value="ttnfitness",tooltip="If X is basename, then X_genes.dat and X_sites.dat will be generated as output files.")
        # can get annotation from GUI; must ask for genome file
        self.value_getters.normalization = create_normalization_input(self.panel, main_sizer) # TTR 
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
            combined_wig = universal.session_data.combined_wigs[0]
            Analysis.inputs.combined_wig = combined_wig.main_path
            
            # 
            # get annotation
            # 

            Analysis.inputs.annotation_path = universal.session_data.annotation # not needed for tnseq_stats
            # FIXME: enable this once I get a valid annotation file example
            # if not transit_tools.validate_annotation(Analysis.inputs.annotation):
            #     return None
            
            # 
            # setup custom inputs
            # 
            for each_key, each_getter in Analysis.instance.value_getters.items():
                try:
                    Analysis.inputs[each_key] = each_getter()
                except Exception as error:
                    raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
            ###transit_tools.log("included_conditions", Analysis.inputs.included_conditions)

            # 
            # get filename for results
            # ttnfitness has 2 output files: "ttnfitness_genes.dat" and "ttnfitness_sites.dat"
            # 

            Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
                default_file_name="ttnfitness_genes.dat",
                output_extensions=u'Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
            )

            Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
                default_file_name="ttnfitness_sites.dat",
                output_extensions=u'Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
            )

            if not Analysis.inputs.output_path: return None
            return Analysis.instance

    @classmethod
    def from_args(cls, args, kwargs): # clean_args() was already called in pytransit/__main__.py
      # usage: python3 %s ttnfitness <comma-separated .wig files> <annotation .prot_table> <genome .fna> <gumbel output file> <output1 file> <output2 file>""" # the old way, with multiple wigs as input

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
            transit_tools.log("Starting tnseq_stats analysis")
            start_time = time.time()

            #######################
            # get data

            if self.inputs.combined_wig!=None: 
              transit_tools.log("Getting Data from %s" % self.inputs.combined_wig)
              sites, data, filenames_in_comb_wig = tnseq_tools.read_combined_wig(self.inputs.combined_wig)

              # read metadata
              # now, select the columns in data corresponding to samples that are replicates of desired condition...

            elif self.inputs.wig_files!=None:
              transit_tools.log("Getting Data")
              (data, position) = transit_tools.get_validated_data( self.inputs.wig_files, wxobj=self.wxobj )

            else: print("error: must provide either combined_wig or list of wig files"); sys.exit(0) ##### use transit.error()?
                
            (K, N) = data.shape

            # normalize the counts
            if self.inputs.normalization and self.inputs.normalization != "nonorm":
              transit_tools.log("Normalizing using: %s" % self.inputs.normalization)
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
    
            transit_tools.log("Getting Genome")
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

            transit_tools.log("processing data")

            TA_sites_df,Models_df,gene_obj_dict = self.calc_ttnfitness(genome,G,self.inputs.gumbel_results_path)

            ###########################
            # write output
            # 
            # note: first header line is filetype, last header line is column headers

            file = sys.stdout # print to console if not output file defined
            if self.inputs.output_path != None:
               file = open(self.inputs.output_path, "w")
            file.write("%s\n" % self.identifier)
            #file.write("#normalization: %s\n" % self.inputs.normalization)

            # I think the column headers will be printed by write_ttnfitness_results(); will it also open the 2 output files?

            self.write_ttnfitness_results(TA_sites_df,Models_df,gene_obj_dict)

            if self.inputs.output_path != None: file.close()
            if universal.interface=="gui" and self.inputs.output_path!=None:
              transit_tools.log(f"Adding File: {self.inputs.output_path}")
              results_area.add(self.inputs.output_path)
            transit_tools.log("Finished TnseqStats")
            transit_tools.log("Time: %0.1fs\n" % (time.time() - start_time))

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

        transit_tools.log("Making Fitness Estimations")
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

        saturation = len(TA_sites_df[TA_sites_df["Insertion Count"] > 0]) / len(
            TA_sites_df
        )
        phi = 1.0 - saturation
        significant_n = math.log10(0.05) / math.log10(phi)

        transit_tools.log("\t + Filtering ES/ESB Genes")
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

        transit_tools.log("\t + Filtering Short Genes. Labeling as Uncertain")
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
        # STLM Predictions
        # transit_tools.log("\t + Making TTN based predictions using loaded STLM")
        # X = filtered_ttn_data.drop(["Orf","Name", "Coord","State","Insertion Count","Local Average","Actual LFC","Upseq TTN","Downseq TTN"],axis=1)
        # X = sm.add_constant(X)
        # model_LFC_predictions = self.STLM_reg.predict(X)
        # filtered_ttn_data["STLM Predicted LFC"]=model_LFC_predictions
        # filtered_ttn_data["STLM Predicted Counts"] = filtered_ttn_data["Local Average"].mul(numpy.power(2,filtered_ttn_data["STLM Predicted LFC"]))

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
        # stlm_predicted_log_counts = numpy.log10(filtered_ttn_data["STLM Predicted Counts"]+0.5)

        Y = numpy.log10(filtered_ttn_data["Insertion Count"] + 0.5)

        transit_tools.log("\t + Fitting M1")
        X1 = pandas.concat([gene_one_hot_encoded, ttn_vectors], axis=1)
        X1 = sm.add_constant(X1)
        results1 = sm.OLS(Y, X1).fit()
        filtered_ttn_data["M1 Pred log Count"] = results1.predict(X1)
        filtered_ttn_data["M1 Predicted Count"] = numpy.power(
            10, (filtered_ttn_data["M1 Pred log Count"] - 0.5)
        )

        # transit_tools.log("\t + Fitting new mod TTN-Fitness")
        # X2 = pandas.concat([gene_one_hot_encoded,stlm_predicted_log_counts],axis=1)
        # X2 = sm.add_constant(X2)
        # results2 = sm.OLS(Y,X2).fit()
        # filtered_ttn_data["mod ttn Pred log Count"] = results2.predict(X2)
        # filtered_ttn_data["mod ttn Predicted Count"] = numpy.power(10, (filtered_ttn_data["mod ttn Pred log Count"]-0.5))

        transit_tools.log("\t + Assessing Models")
        # create Models Summary df
        Models_df = pandas.DataFrame(results1.params[1:-256], columns=["M1 Coef"])
        Models_df["M1 Pval"] = results1.pvalues[1:-256]
        Models_df["M1 Adjusted Pval"] = statsmodels.stats.multitest.fdrcorrection(
            results1.pvalues[1:-256], alpha=0.05
        )[1]
        # Models_df["mod ttn Coef"] = results2.params[1:-1]
        # Models_df["mod ttn Pval"] = results2.pvalues[1:-1]
        # Models_df["mod ttn Adjusted Pval"] = statsmodels.stats.multitest.fdrcorrection(results2.pvalues[1:-1],alpha=0.05)[1]

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

        # mask using mod TTN fitness
        # Models_df.loc[(Models_df["mod ttn Coef"]>0) & (Models_df["mod ttn Adjusted Pval"]<0.05),"mod ttn States"]="GA"
        # Models_df.loc[(Models_df["mod ttn Coef"]<0) & (Models_df["mod ttn Adjusted Pval"]<0.05),"mod ttn States"]="GD"
        # Models_df.loc[(Models_df["mod ttn Coef"]==0) & (Models_df["mod ttn Adjusted Pval"]<0.05),"mod ttn States"]="NE"
        # Models_df.loc[(Models_df["mod ttn Adjusted Pval"]>0.05),"mod ttn States"]="NE"

        return (TA_sites_df,Models_df,gene_obj_dict)

    def write_ttnfitness_results(TA_sites_df,Models_df,gene_obj_dict):
        transit_tools.log("Writing To Output Files")
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
                # TA_sites_df[TA_sites_df["Coord"].isin(coords_orf)]['mod ttn Predicted Count'] = filtered_ttn_data[filtered_ttn_data["Coord"].isin(coords_orf)]["mod ttn Predicted Count"]
            # M1 info
            if "_" + g in Models_df.index:
                M1_coef = Models_df.loc["_" + g, "M1 Coef"]
                M1_adj_pval = Models_df.loc["_" + g, "M1 Adjusted Pval"]
                modified_M1 = math.exp(
                    M1_coef - statistics.median(Models_df["M1 Coef"].values.tolist())
                )
                # mod_M1_coef = Models_df.loc["_"+g,"mod ttn Coef"]
                # mod_M1_adj_pval = Models_df.loc["_"+g,"mod ttn Adjusted Pval"]
                # mod_modified_M1 = math.exp(mod_M1_coef - statistics.median(Models_df["mod ttn Coef"].values.tolist()))
            else:
                M1_coef = None
                M1_adj_pval = None
                modified_M1 = None
                # mod_M1_coef = None
                # mod_M1_adj_pval = None
                # mod_modified_M1 = None

            # States
            gumbel_bernoulli_call = gumbel_bernoulli_gene_calls[g]
            if gumbel_bernoulli_call == "E":
                gene_ttn_call = "ES"
                # mod_gene_ttn_call = "ES"
            elif gumbel_bernoulli_call == "EB":
                gene_ttn_call = "ESB"
                # mod_gene_ttn_call = "ESB"
            else:
                if "_" + g in Models_df.index:
                    gene_ttn_call = Models_df.loc["_" + g, "Gene+TTN States"]
                    # mod_gene_ttn_call = Models_df.loc["_"+g,"mod ttn States"]
                else:
                    gene_ttn_call = "U"  # these genes are in the uncertain genes list
                    # mod_gene_ttn_call = "U"
            TA_sites_df.loc[
                (TA_sites_df["Orf"] == g), "TTN-Fitness Assessment"
            ] = gene_ttn_call
            # TA_sites_df.loc[(TA_sites_df["Orf"]==g), 'Mod TTN-Fitness Assessment'] = mod_gene_ttn_call
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
        # mod_assesment_cnt = output_df["Mod TTN-Fitness Assessment"].value_counts()

        self.output.write("#TTNFitness\n")
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
                    ",".join(self.ctrldata).encode("utf-8"),
                    self.annotation_path.encode("utf-8"),
                    self.output.name.encode("utf-8"),
                )
            )
        else:
            self.output.write("#Console: python3 %s\n" % " ".join(sys.argv))

        self.output.write("#Data: %s\n" % (",".join(self.ctrldata).encode("utf-8")))
        self.output.write(
            "#Annotation path: %s\n" % self.annotation_path.encode("utf-8")
        )
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#Saturation of Dataset: %s\n" % (saturation))
        self.output.write(
            "#Assesment Counts: %s ES, %s ESB, %s GD, %s GA, %s NE, %s U \n"
            % (
                assesment_cnt["ES"],
                assesment_cnt["ESB"],
                assesment_cnt["GD"],
                assesment_cnt["GA"],
                assesment_cnt["NE"],
                assesment_cnt["U"],
            )
        )
        # self.output.write("#Mod Assesment Counts: %s ES, %s ESB, %s GD, %s GA, %s NE, %s U \n" % (mod_assesment_cnt["ES"],mod_assesment_cnt["ESB"],mod_assesment_cnt["GD"],mod_assesment_cnt["GA"],mod_assesment_cnt["NE"],mod_assesment_cnt["U"]))

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

        output2_data = TA_sites_df.to_csv(header=True, sep="\t", index=False).split(
            "\n"
        )
        vals = "\n".join(output2_data)
        self.output2_file.write(vals)
        self.output2_file.close()

        output_data = output_df.to_csv(header=True, sep="\t", index=False).split("\n")
        vals = "\n".join(output_data)
        self.output.write(vals)
        self.output.close()

        transit_tools.log("")  # Printing empty line to flush stdout
        transit_tools.log("Adding File: %s" % (self.output.name))
        results_area.add(self.output.name)
        self.finish()
        transit_tools.log("Finished TTNFitness Method")


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
            name=transit_tools.basename(self.path),
            type=Analysis.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(title="tnseq_stats",heading="",column_names=self.column_names,rows=self.rows).Show(),
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
