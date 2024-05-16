import sys
import os
import time
import ntpath
import math
import random
import datetime
import collections
import heapq
import numpy as np
import gzip

import numpy

from pytransit.generic_tools import csv, misc, informative_iterator
from pytransit.specific_tools import gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.specific_tools.transit_tools import wx, basename
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.generic_tools.informative_iterator import ProgressBar

@misc.singleton
class Method:
    name = "CGI"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name} - use CRISPRi-DR to analyze Chemical-Genetic Interactions"
    description = f"""Perform {name} analysis"""
    
    valid_cli_flags = [
        "-use_negatives", #use negative controls to calculate pvalues
        "--fixed",  # fixed axes
        "-origx", # non-log axes
        "-origy", # non-log axes
        "-no_uninduced", #no uninduced ATC values
        "-delete_temp_fastQ", #if gz files, this flag indicates whether user would like to delete the temp files
    ]
    usage_string = f"""
        Usage (6 sub-commands):

    {console_tools.subcommand_prefix} {cli_name} extract_counts <fastq file> <sgRNA info file> <barcode_col> <output counts file>
        Optional Arguements:
            -delete_temp_fastQ := if fast files are gz files, this flag indicates whether user would like to delete the temp files
    
    {console_tools.subcommand_prefix} {cli_name} create_combined_counts <comma seperated headers> <counts file 1> <counts file 2> ... <counts file n> <combined counts file>
    
    {console_tools.subcommand_prefix} {cli_name} extract_abund <combined counts file> <metadata file> <control condition> <sgRNA info file> <efficacy column name> <orf column name> <uninduced ATC file> <drug> <days>  <fractional abundance file>
        Optional Arguments:
            -no_uninduced := flag to calculate fractional abundances without input uninduced ATC file. If this flag is set, then the uninduced ATC file should not be provided

    {console_tools.subcommand_prefix} {cli_name} run_model <fractional abundance file>  <CRISPRi DR results file>
        Optional Arguments:
            -use_negatives := flag to use negative controls to calculate significance of coefficients of concentration dependence

    {console_tools.subcommand_prefix} {cli_name} visualize <fractional abundance> <gene> <output figure location> [Optional Arguments]
        Optional Arguments: 
            --fixed xmin=x,xmax=x,ymin=y,ymax=y := set the values you would to be fixed in this comma seperated format. Not all values need to be set for ex, a valid argument is "xmin=0,ymax=5"
            -origx := flag to turn on original scale axes rather than log scale for Concentration default=off
            -origy := flag to turn on original scale axes rather than log scale for Realtive Abundances default=off
    """.replace("\n        ", "\n")
    
    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        logging.log("Please provide a subcommand for cgi")
        print(Method.usage_string)

    @staticmethod
    @cli.add_command(cli_name, "extract_counts")
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=3)

        fastq_file= args[0]
        sgRNA_info_file = args[1]
        barcode_col = args[2]
        counts_file = args[3]
        delete_temp_fastQ = True if "delete_temp_fastQ" in kwargs else False
        Method.extract_counts(fastq_file, sgRNA_info_file, barcode_col, counts_file, delete_temp_fastQ= delete_temp_fastQ)

    @staticmethod
    @cli.add_command(cli_name, "create_combined_counts")
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        #console_tools.enforce_number_of_args(args, Method.usage_string, greater=2)

        headers = args[0].split(",")
        counts_file_list = args[1:-1]
        combined_counts_file = args[-1]
        Method.create_combined_counts(headers,counts_file_list, combined_counts_file)


    @staticmethod
    @cli.add_command(cli_name, "extract_abund")
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=7)

        combined_counts_file      = args[0]
        metadata_file             = args[1]
        control_condition         = args[2]
        sgRNA_info_file      = args[3]
        efficacy_col = args[4]
        orf_col = args[5]
        no_uninduced = True if "no_uninduced" in kwargs else False

        if no_uninduced:
            uninduced_atc_file        = None
            drug                      = args[6]
            days                      = args[7]
            fractional_abundance_file = args[8]
        else:
            uninduced_atc_file        = args[6]
            drug                      = args[7]
            days                      = args[8]
            fractional_abundance_file = args[9]

        Method.extract_abund(
            combined_counts_file,
            metadata_file,
            control_condition,
            sgRNA_info_file,
            efficacy_col,
            orf_col,
            uninduced_atc_file,
            drug,
            days,
            fractional_abundance_file,
            no_uninduced = no_uninduced
        )

    @staticmethod
    @cli.add_command(cli_name, "run_model")
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=2)

        ifile_path = args[0]
        cdr_output_file = args[1]
        use_negatives = True if "use_negatives" in kwargs else False
        Method.run_model(ifile_path, cdr_output_file, use_negatives=use_negatives)

    @staticmethod
    @cli.add_command(cli_name, "visualize")
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=3)

        frac_abund_file= args[0]
        gene = args[1]
        fig_location = args[2]

        fixed_axes = kwargs.get("fixed", "") # did it this way rather than having xmin, xmax etc. as seperate flag because the - in say -2.5 would trigger flag recognition
        origx_axis = True if "origx" in kwargs else False
        origy_axis = True if "origy" in kwargs else False

        Method.visualize(frac_abund_file, gene, fig_location, fixed = fixed_axes, origx = origx_axis, origy=origy_axis)

    def reverse_complement(self, seq):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        s = list(seq)
        s.reverse()
        for i in range(len(s)):
            s[i] = complement.get(s[i], s[i])  # if unknown, leave as it, e.g > or !
        s = "".join(s)
        return s

    def extract_counts(self, fastq_file_initial, sgRNA_info_file, barcode_col, counts_file, delete_temp_fastQ= False):
        import pandas as pd

        is_gz_file = False
        if fastq_file_initial[-3:]==".gz": 
            is_gz_file = True 
            logging.log("%s Detected as a gzipped file \n" % fastq_file_initial)
        
        if is_gz_file:
            logging.log("Uncompressing %s \n" % fastq_file_initial)
            fastq_file = "temp_"+fastq_file_initial[0:-3]
            outfil = open(fastq_file, "wb+")
            for line in gzip.open(fastq_file_initial):
                outfil.write(line)
        else:
            fastq_file = fastq_file_initial


        f = open(counts_file, "w")
        IDs = []
        barcodemap = {}  # hash from barcode to full ids
        sgRNA_info_df = pd.read_csv(sgRNA_info_file, sep="\t", index_col=0)
        for id, row in sgRNA_info_df.iterrows():
            barcode= row[barcode_col]
            IDs.append(id)
            barcodemap[self.reverse_complement(barcode)] = id

        counts = {}
        A, B = "GTACAAAAAC", "TCCCAGATTA"
        lenA = len(A)
        cnt, nreads, recognized = 0, 0, 0
        for line in open(fastq_file):
            cnt += 1
            if cnt % 4 == 2:
                nreads += 1
                if nreads % 1000000 == 0:
                    logging.log(
                        "reads=%s, recognized barcodes=%s (%0.1f%%)\n"
                        % (nreads, recognized, 100.0 * recognized / float(nreads))
                    )
                seq = line.rstrip()
                a = seq.find(A)
                if a == -1:
                    continue
                b = seq.find(B)
                if b == -1:
                    continue
                sz = b - (a + lenA)
                if sz < 10 or sz > 30:
                    continue
                barcode = seq[
                    a + lenA : b
                ]  # these are reverse-complements, but rc(barcodes) stored in hash too
                if barcode not in barcodemap:
                    continue
                id = barcodemap[barcode]
                if id not in counts:
                    counts[id] = 0
                counts[id] += 1
                recognized += 1

        for id in IDs:
            vals = [id, counts.get(id, 0)]
            f.write("\t".join([str(x) for x in vals]))
            f.write("\n")
        f.close()

        if is_gz_file and delete_temp_fastQ:
            os.remove(fastq_file)

    def create_combined_counts(self, headers, counts_list, combined_counts_file):
        import pandas as pd

        df_list = []
        for f in counts_list:
            logging.log("Adding in file # %s \n" % f)
            counts_df = pd.read_csv(f, sep="\t", comment="#", index_col=0)
            df_list.append(counts_df)
        combined_df = pd.concat(df_list, axis=1)
        combined_df.columns = headers
        combined_df.to_csv(combined_counts_file, sep="\t")
        logging.log(
            "Number of sgRNAs in combined counts file (present in all counts files): %d \n"
            % len(combined_df)
        )

    def extract_abund(
        self,
        combined_counts_file,
        metadata_file,
        control_condition,
        sgRNA_info_file,
        efficacy_col,
        orf_col,
        uninduced_atc_file,
        drug,
        days,
        fractional_abundance_file,
        PC=1e-8,
        no_uninduced=False,
    ):
        import pandas as pd

        metadata = pd.read_csv(metadata_file, sep="\t", comment="#")
        available_drugs = misc.no_duplicates(metadata["drug"].values.tolist())
        available_days = misc.no_duplicates(
            metadata["days_predepletion"].values.tolist()
        )

        # Misnamed args check
        assert drug in available_drugs, logging.error(
            f"{repr(drug)} is not one of the valid drug options from your metadata file ({repr(available_drugs)}). Add the drug's information in the metadata file or select a different drug"
        )
        assert int(days) in available_days, logging.error(
            f"{repr(days)} is not one of the valid 'days' options in your metadata days of predepletion column ({repr(available_days)}). Add the day's information in the metadata file or select a different day"
        )
        assert control_condition in available_drugs, logging.error(
            f"{repr(control_condition)} is not one of the valid control_condition options from your metadata file ({repr(available_drugs)}). Add the corresponding information in the metadata file or select a different control"
        )

        metadata = metadata[
            ((metadata["drug"] == drug) | (metadata["drug"] == control_condition))
            & (metadata["days_predepletion"] == int(days))
        ]
        available_drugs = misc.no_duplicates(metadata["drug"].values.tolist())
        available_days = misc.no_duplicates(
            metadata["days_predepletion"].values.tolist()
        )

        # Bad combinations check
        assert len(metadata) != 0, logging.error(
            f"This combination (drug={drug}, control_condition={control_condition}, days={days}) of conditions does not exist in your metadata file. Please select one that does"
        )
        assert drug in available_drugs, logging.error(
            f"{repr(drug)} is not one of the valid options from your metadata file ({repr(available_drugs)}). Add the drug's information in the metadata file or select a different drug"
        )
        assert int(days) in available_days, logging.error(
            f"{repr(days)} is not found in your metadata days of predepletion column (repr({available_days})). Add the day's information in the metadata file or select a different day"
        )
        assert control_condition in available_drugs, logging.error(
            f"{repr(control_condition)} is not one of the valid options from your metadata file ({repr(available_drugs)}). Add the corresponding information in the metadata file or select a different control"
        )

        metadata = metadata.sort_values(by=["conc_xMIC"])
        column_names = metadata["column_name"].values.tolist()
        concs_list = metadata["conc_xMIC"].values.tolist()
        
        logging.log("# Condition Tested : "+str(drug)+"-"+str(days))
        headers = []
        combined_counts_df = pd.read_csv(combined_counts_file,sep="\t", index_col=0)
        combined_counts_df = combined_counts_df[column_names]
        combined_counts_df = combined_counts_df[~combined_counts_df.index.str.contains("Empty")]


        if(len(combined_counts_df.columns)==0):
            logging.error("The samples associated with the selected drugs do not exist in your combined counts file. Please select one that does and check your metadata file has corresponding column names")
        elif(len(combined_counts_df.columns)<len(metadata)):
            logging.log("WARNING: Not all of the samples from the metadata based on this criteron have a column in the combined counts file")
      
        sgrna_info_df = pd.read_csv(sgRNA_info_file,sep="\t", index_col=0)
        if efficacy_col=="":
            sgrna_efficacy = sgrna_info_df[sgrna_info_df.columns[2]]
        else:
            sgrna_efficacy = sgrna_info_df[[efficacy_col]]
        sgrna_efficacy.columns = ["sgRNA efficacy"]

        if orf_col=="":
            orf_info = sgrna_info_df[sgrna_info_df.columns[0]]
        else:  
            orf_info = sgrna_info_df[[orf_col]]
        orf_info.columns = ["Orf"]

        if no_uninduced:
            PC = 1e-8
            logging.log("Calculating Fractional Abundances without Uninduced Abundances")
            conc0_cols = metadata[metadata["conc_xMIC"]==0]["column_name"].values.tolist()
            no_dep_df = pd.DataFrame(index=combined_counts_df.index)
            no_dep_df["SCV"] = (combined_counts_df.std(axis=1)+PC)/(combined_counts_df.mean(axis=1)+PC)
            no_dep_df["Conc0 Mean"] = combined_counts_df[conc0_cols].mean(axis=1)
            no_dep_df["uninduced ATC values"] = no_dep_df["Conc0 Mean"] * np.exp(2*no_dep_df["SCV"])
            no_dep_df = no_dep_df[["uninduced ATC values"]]
            no_dep_df["uninduced ATC values"] = (no_dep_df["uninduced ATC values"] / no_dep_df["uninduced ATC values"].sum())
        else:
            logging.log("Calculating Fractional Abundances with Uninduced Abundances")
            no_dep_df = pd.read_csv(uninduced_atc_file, sep="\t", index_col=0, header=None, comment="#")
            no_dep_df = no_dep_df.iloc[:, -1:]
            no_dep_df.columns = ["uninduced ATC values"]
            no_dep_df["uninduced ATC values"] = (no_dep_df["uninduced ATC values"] / no_dep_df["uninduced ATC values"].sum())

        abund_df = pd.concat([orf_info, sgrna_efficacy, no_dep_df,combined_counts_df], axis=1)
        abund_df= abund_df[~abund_df.index.str.contains("Empty")]
        logging.log("Disregarding Empty sgRNAs")
        logging.log("%d usable sgRNAs (including Negatives) are all of the following files : sgRNA efficacy metadata, uninduced ATC counts file, combined counts file"%len(abund_df))
        abund_df["uninduced ATC values"] = abund_df["uninduced ATC values"].replace(np.nan,1)


        headers = ["Orf","sgRNA efficacy", "uninduced ATC values"]
        if fractional_abundance_file != None:
            f = open(fractional_abundance_file, "w")
        for i, col in enumerate(column_names):
            abund_df[col] = abund_df[col] / abund_df[col].sum()
            abund_df[col] = (abund_df[col] + PC) / (abund_df["uninduced ATC values"] + PC)
            headers.append(str(concs_list[i]) + "_" + str(i))
            if fractional_abundance_file != None:
                f.write("# " + str(concs_list[i]) + " conc_xMIC" + " - " + col + "\n")

        abund_df.columns = headers
        
        #Lines to ensure first column name (index column) has sgRNA label and any orf names will nulls are dropped
        abund_df = abund_df[~abund_df["Orf"].isna()]
        abund_df["sgRNA"] = abund_df.index
        abund_df.set_index("sgRNA", inplace=True)
        
        if fractional_abundance_file != None:
            abund_df.to_csv(f, sep="\t")
        
        return abund_df



    def run_model(self, frac_abund_file, cdr_output_file, use_negatives=False):
        import pandas as pd
        import numpy as np
        from mne.stats import fdr_correction
        import statsmodels.api as sm
        from pytransit.components import parameter_panel

        transit_tools.require_r_to_be_installed(required_r_packages=[ "locfdr" ])
        
        frac_abund_df = pd.read_csv(frac_abund_file, sep="\t",comment='#')
        if use_negatives:
            logging.log("Alert : -use_negatives flag has been provided, significance of genes will be calculated using Negative controls")
            frac_abund_df= frac_abund_df[~(frac_abund_df["sgRNA"].str.contains("Empty"))]
            frac_abund_df = frac_abund_df.fillna(1)
            total_neg_genes= int(np.sqrt(len(frac_abund_df[frac_abund_df["Orf"].str.contains("Negative")])))
            counter =0
            for neg_sgRNA in set(frac_abund_df[frac_abund_df["Orf"].str.contains("Negative")]["sgRNA"]):
                neg_gene = counter%total_neg_genes
                counter= counter+1
                frac_abund_df.loc[frac_abund_df["sgRNA"]==neg_sgRNA,"Orf"] = "Negative_"+str(neg_gene)
        else: 
            frac_abund_df= frac_abund_df[~(frac_abund_df["sgRNA"].str.contains("Negative") | frac_abund_df["sgRNA"].str.contains("Empty"))]
        
        frac_abund_df = frac_abund_df.dropna()
        drug_output = []
        for i,orf in enumerate(set(frac_abund_df["Orf"])):
            logging.log("Analyzing Gene #"+str(i)+" - "+str(orf))
            gene_df = frac_abund_df[frac_abund_df["Orf"]==orf]
            gene_df = gene_df.drop(columns=["Orf", "uninduced ATC values"])

            melted_df = gene_df.melt(
                id_vars=["sgRNA", "sgRNA efficacy"], var_name="conc", value_name="abund"
            )
            try:
                melted_df["conc"] = (
                    melted_df["conc"].str.split("_", expand=True)[0].astype(float)
                )
                min_conc = min(melted_df[melted_df["conc"] > 0]["conc"])
            except Exception as error:
                import code; code.interact(local={**globals(),**locals()})
            melted_df.loc[melted_df["conc"] == 0, "conc"] = min_conc / 2
            melted_df["abund"] = [
                0.01
                + (1 - 0.01) * (1 - np.exp(-2 * float(i))) / (1 + np.exp(-2 * float(i)))
                for i in melted_df["abund"]
            ]
            melted_df["logsig abund"] = [
                np.nan if (1 - x) == 0 else np.log10(float(x) / (1 - float(x)))
                for x in melted_df["abund"]
            ]
            melted_df["log conc"] = [np.log2(float(x)) for x in melted_df["conc"]]

            melted_df = melted_df.dropna()
            if len(melted_df.index) < 2:
                print("Found a small gene", orf)
                drug_output.append([orf, len(gene_df)] + [np.nan] * 6)
                continue
            status = " "
            if len(gene_df)<5:
                logging.log("NOTE: ",orf," has less than 5 sgRNAs")
                status = "< 5 sgRNAs"

            Y = melted_df["logsig abund"]
            X = melted_df.drop(columns=["abund", "logsig abund", "sgRNA", "conc"])
            X = sm.add_constant(X)
            model = sm.OLS(Y, X)
            results = model.fit()
            coeffs = results.params
            pvals = results.pvalues
            drug_output.append(
                [status, orf, len(gene_df)]
                + coeffs.values.tolist()
                + pvals.values.tolist()
            )
            percentage = (100.0 * i / len(set(frac_abund_df["Orf"])))
            if gui.is_active:
                parameter_panel.progress_update(f"Analyzing Genes  {percentage:.2f}%\r", percentage)

        drug_out_df = pd.DataFrame(
            drug_output,
            columns=[
                "Status",
                "Orf",
                "Nobs",
                "intercept",
                "coefficient sgrna efficacy",
                "coefficient concentration dependence",
                "pval intercept",
                "pval sgrna efficacy",
                "pval concentration dependence",
            ],
        )

        drug_out_df["intercept"] = round(drug_out_df["intercept"],6)
        drug_out_df["coefficient sgrna efficacy"] = round(drug_out_df["coefficient sgrna efficacy"],6)
        drug_out_df["coefficient concentration dependence"] = round(drug_out_df["coefficient concentration dependence"],6)
        drug_out_df["pval intercept"] = round(drug_out_df["pval intercept"],6)
        drug_out_df["pval sgrna efficacy"] = round(drug_out_df["pval sgrna efficacy"],6)
        drug_out_df["pval concentration dependence"] = round(drug_out_df["pval concentration dependence"],6)
    
        mask = np.isfinite(drug_out_df["pval concentration dependence"])
        pval_corrected = np.full(
            drug_out_df["pval concentration dependence"].shape, np.nan
        )
        pval_corrected[mask] = fdr_correction(
            drug_out_df["pval concentration dependence"][mask]
        )[1]
        drug_out_df["qval concentration dependence"] = pval_corrected
        drug_out_df = drug_out_df.replace(np.nan, 1)



        if use_negatives==True:
            negatives_output = drug_out_df[drug_out_df["Orf"].str.contains("Negative")]
            neg_mean = negatives_output["coefficient concentration dependence"].mean()
            neg_stdev = negatives_output["coefficient concentration dependence"].std()
            drug_out_df["Z score of concentration dependence"] = (drug_out_df["coefficient concentration dependence"] - neg_mean)/neg_stdev
            drug_out_df["qval concentration dependence"] = round(drug_out_df["qval concentration dependence"],6)
            drug_out_df["Z score of concentration dependence"] = round(drug_out_df["Z score of concentration dependence"],6)

            drug_out_df["Significant Interactions"] = [0] * len(drug_out_df)
            drug_out_df.loc[(drug_out_df["qval concentration dependence"]<0.05) & (drug_out_df["Z score of concentration dependence"]<-2),"Significant Interactions"]=-1
            drug_out_df.loc[(drug_out_df["qval concentration dependence"]<0.05) &  (drug_out_df["Z score of concentration dependence"]>2),"Significant Interactions"]=1

            drug_out_df.insert(0, "Significant Interactions", drug_out_df.pop("Significant Interactions"))

            n = len(drug_out_df[drug_out_df["Significant Interactions"] != 0])
            depl_n = len(drug_out_df[drug_out_df["Significant Interactions"] == -1])
            enrich_n = len(drug_out_df[drug_out_df["Significant Interactions"] == 1])
            logging.log("%d Total Significant Gene Interactions" % n)
            logging.log("%d Significant Gene Depletions" % depl_n)
            logging.log("%d Significant Gene Enrichments" % enrich_n)


        else:
            drug_out_df["Z score of concentration dependence"] = (drug_out_df["coefficient concentration dependence"] - drug_out_df["coefficient concentration dependence"].mean())/drug_out_df["coefficient concentration dependence"].std()
        
            drug_out_df["qval concentration dependence"] = round(drug_out_df["qval concentration dependence"],6)
            drug_out_df["Z score of concentration dependence"] = round(drug_out_df["Z score of concentration dependence"],6)

            ## Empirical Bayes Adjustment 
            drug_out_df.to_csv("./temp.txt", sep="\t", index=False)

        
            from pytransit.specific_tools.transit_tools import r, globalenv
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri

            pandas2ri.activate()

            r('''
            library('locfdr')
            locfdr_R <- function(verbose=TRUE) {
                data = read.table("./temp.txt", sep='\t',head=T)

                lim=10
                z = data$Z.score.of.concentration.dependence
                z[z > lim] = lim
                z[z < -lim] = -lim
            
                mod = locfdr(z,nulltype=1, df=20)
                data$locfdr = mod$fdr

                temp = data
                temp = temp[order(temp$locfdr),]

                EFDR = NULL
                for (i in 1:nrow(temp))
                {
                efdr = mean(temp$locfdr[1:i])
                EFDR = rbind(EFDR,efdr)
                }
                temp$ebFDR = EFDR
                data$ebFDR = temp[rownames(data),"ebFDR"]
            
                write.table(data,"./temp.txt",sep="\t",row.names=F,quote=F)
            }
            ''')
            r_f = globalenv['locfdr_R']
            r_f()
            drug_out_df = pd.read_csv("./temp.txt", sep="\t")
            drug_out_df.columns = [c.replace(".", " ") for c in drug_out_df.columns]
            

            import os
            os.remove("./temp.txt")
        #################

            drug_out_df["Significant Interactions"] = [0] * len(drug_out_df)
            drug_out_df.loc[(drug_out_df["qval concentration dependence"]<0.05) & (drug_out_df["ebFDR"]<0.05)& (drug_out_df["coefficient concentration dependence"]<0),"Significant Interactions"]=-1
            drug_out_df.loc[(drug_out_df["qval concentration dependence"]<0.05) & (drug_out_df["ebFDR"]<0.05)& (drug_out_df["coefficient concentration dependence"]>0),"Significant Interactions"]=1

            drug_out_df.insert(0, "Significant Interactions", drug_out_df.pop("Significant Interactions"))

            n = len(drug_out_df[drug_out_df["Significant Interactions"] != 0])
            depl_n = len(drug_out_df[drug_out_df["Significant Interactions"] == -1])
            enrich_n = len(drug_out_df[drug_out_df["Significant Interactions"] == 1])
            logging.log("%d Total Significant Gene Interactions" % n)
            logging.log("%d Significant Gene Depletions" % depl_n)
            logging.log("%d Significant Gene Enrichments" % enrich_n)

        #drug_out_df = drug_out_df.replace(r"\s+", np.nan, regex=True).replace("", np.nan)
        try:
            transit_tools.write_result(
                path=cdr_output_file,  # path=None means write to STDOUT
                file_kind=Method.identifier,
                rows=tuple(each.tolist() for each in drug_out_df.iloc),
                column_names=drug_out_df.columns,
                extra_info=dict(
                    stats=dict(),
                    parameters=dict(
                        frac_abund_file=frac_abund_file,
                        # TODO: add the other parameters
                    ),
                    results=dict(
                        total_hits = n,
                        depleted_hits = depl_n,
                        enriched_hits = enrich_n
                    )
                ),
            )
        except Exception as error:
            print("error", error)
            import code; code.interact(local={**globals(),**locals()})



    def visualize(
        self, fractional_abundances_file, gene, fig_location, fixed, origx, origy
    ):
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import numpy as np
        import statsmodels.api as sm
        import re
        
        # if given multiple genes, call self with each
        if "," in gene:
            for each_gene in gene.split(","):
                each_gene = each_gene.strip()
                self.visualize(
                    fractional_abundances_file=fractional_abundances_file,
                    gene=each_gene,
                    fig_location=fig_location,
                    fixed=fixed,
                    origx=origx,
                    origy=origy,
                )
            return

        abund_df = pd.read_csv(fractional_abundances_file, sep="\t", comment="#")
        
        abund_df_lowercased = pd.DataFrame(abund_df)
        #abund_df_lowercased["gene"] = tuple(str(each).lower() for each in abund_df["gene"].values.tolist())
        abund_df_lowercased["Orf"] = tuple(str(each).lower() for each in abund_df["Orf"].values.tolist())
        gene_no_case = gene.lower()
        abund_df_filtered = abund_df_lowercased[(abund_df_lowercased["Orf"] == gene_no_case)]
        if len(abund_df_filtered) == 0:
            available_genes = list(abund_df["Orf"].tolist())
            logging.error(f"Gene {gene}, was not one of the available genes: {repr(available_genes)}")
        abund_df = abund_df_filtered.reset_index(drop=True)
        all_slopes = []
        
        

        df_list = []
        for idx, row in abund_df.iterrows():
            try:
                logging.log("Fitting sgRNA # : %d \n" % idx)
                raw_Y = row[5:].values
                Y = [max(0.01, x) for x in raw_Y]
                Y = [np.log10(x) for x in Y]
            
                raw_X = abund_df.columns[5:]
                raw_X = [float(i.split("_")[0]) for i in raw_X]
                min_conc = min([i for i in raw_X if i > 0])
                X = [min_conc / 2 if i == 0 else i for i in raw_X]
                X = [np.log2(float(x)) for x in X]

                data = pd.DataFrame(
                    {"Log (Concentration)": X, "Log (Relative Abundance)": Y}
                )
                X = pd.DataFrame({"log concentration": X})
                X_in = sm.add_constant(X, has_constant="add")
                results = sm.OLS(Y, X_in).fit()
                all_slopes.append(results.params[1])
                data["sgRNA efficacy"] = [row["sgRNA efficacy"]] * len(data)
                data["slope"] = [results.params[1]] * len(data)
                data["Concentration"] = raw_X
                data["Relative Abundance"] = raw_Y
                df_list.append(data)
            except Exception as error:
                print(error)
                import code; code.interact(local={**globals(),**locals()})

        plot_df = pd.concat(df_list)

        cmap = sns.color_palette("Spectral", as_cmap=True)
        #palette = [mpl.colors.rgb2hex(cmap(i)) for i in range(len(abund_df))]
        palette  = sns.color_palette("Spectral", len(abund_df))

        if origx == True:
            xmin = plot_df["Concentration"].min()
            xmax = plot_df["Concentration"].max()
        else:
            xmin = plot_df["Log (Concentration)"].min()
            xmax = plot_df["Log (Concentration)"].max()

        if origy == True:
            ymin = plot_df["Relative Abundance"].min()
            ymax = plot_df["Relative Abundance"].max()
        else:
            ymin = plot_df["Log (Relative Abundance)"].min()
            ymax = plot_df["Log (Relative Abundance)"].max()

        if origx == False and origy == False:
            g = sns.lmplot(
                data=plot_df,
                x="Log (Concentration)",
                y="Log (Relative Abundance)",
                hue="sgRNA efficacy",
                palette=palette,
                legend=False,
                ci=None,
                scatter=False,
                line_kws={"lw": 0.75},
            )
        elif origx == False and origy == True:
            g = sns.lmplot(
                data=plot_df,
                x="Log (Concentration)",
                y="Relative Abundance",
                hue="sgRNA efficacy",
                palette=palette,
                legend=False,
                ci=None,
                scatter=False,
                line_kws={"lw": 0.75},
            )
        elif origx == True and origy == False:
            g = sns.lmplot(
                data=plot_df,
                x="Concentration",
                y="Log (Relative Abundance)",
                hue="sgRNA efficacy",
                palette=palette,
                legend=False,
                ci=None,
                scatter=False,
                line_kws={"lw": 0.75},
            )
        else:
            g = sns.lmplot(
                data=plot_df,
                x="Concentration",
                y="Relative Abundance",
                hue="sgRNA efficacy",
                palette=palette,
                legend=False,
                ci=None,
                scatter=False,
                line_kws={"lw": 0.75},
            )

        range_values = fixed.split(",")

        xmin_list = [x for x in range_values if re.search("xmin=", x)]
        if len(xmin_list) > 1:
            logging.log(
                "You have provided more than one xmin value. Figure will be created using flexible xmin"
            )
        elif len(xmin_list) == 1:
            xmin_temp = float(xmin_list[0].split("=")[1])
            if xmin_temp > xmax:
                logging.log(
                    "You have provided an xmin value greater than the xmax value. Figure will be created using flexible xmin"
                )
            else:
                xmin = xmin_temp

        xmax_list = [x for x in range_values if re.search("xmax=", x)]
        if len(xmax_list) > 1:
            logging.log(
                "You have provided more than one xmax value. Figure will be created using flexible xmax"
            )
        elif len(xmax_list) == 1:
            xmax_temp = float(xmax_list[0].split("=")[1])
            if xmax_temp < xmin:
                logging.log(
                    "You have provided an xmax value less than the xmin value. Figure will be created using flexible xmax"
                )
            else:
                xmax = xmax_temp

        ymin_list = [x for x in range_values if re.search("ymin=", x)]
        if len(ymin_list) > 1:
            logging.log(
                "You have provided more than one ymin value. Figure will be created using flexible ymin"
            )
        elif len(ymin_list) == 1:
            ymin_temp = float(ymin_list[0].split("=")[1])
            if ymin_temp < ymax:
                logging.log(
                    "You have provided an ymin value greater than the ymax value. Figure will be created using flexible ymin"
                )
            else:
                ymin = ymin_temp

        ymax_list = [x for x in range_values if re.search("ymax=", x)]
        if len(ymax_list) > 1:
            logging.log(
                "You have provided more than one ymax value. Figure will be created using flexible ymax"
            )
        elif len(ymax_list) == 1:
            ymax_temp = float(ymax_list[0].split("=")[1])
            if ymax_temp < ymin:
                logging.log(
                    "You have provided an ymax value less than the ymin value. Figure will be created using flexible ymax"
                )
            else:
                ymax = ymax_temp

        sm1 = mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(
                vmin=abund_df["sgRNA efficacy"].min(), vmax=abund_df["sgRNA efficacy"].max(), clip=False
            ),
            cmap=cmap,
        )
        g.figure.colorbar(sm1, shrink=0.8, aspect=50, label="sgRNA efficacy")
        g.set(ylim=(ymin, ymax), xlim=(xmin, xmax))
        plt.gca().set_title(gene, wrap=True)
        plt.tight_layout()
        if not fig_location:
            plt.show()
        else:
            plt.savefig(fig_location+".png")

    @gui.add_menu("Method", "CRISPRi", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)

    def define_panel(self, _):
        from pytransit.components import panel_helpers, parameter_panel
        import pandas as pd

        with panel_helpers.NewPanel() as (panel, main_sizer):
            parameter_panel.set_instructions(
                title_text=self.name,
                sub_text="",
                method_specific_instructions="""
                    The CRISPRi-DR methods is designed to analyze CRISPRi libraries from Chemical-Genetic Interation (CGI) experiments and identify significant CGIs, i.e. genes that affect sensitivity to the drug when depleted.  Dropdown boxes to select Drug and Control will appear after the metadata is loaded.
                """.replace(
                    "\n                    ", "\n"
                ),
            )
            self.value_getters = LazyDict(
                control_condition=lambda : None,
                drug=lambda : None,
                days=lambda : None,
                cgi_folder=lambda : None,
                combined_counts_file=lambda : None,
                metadata_file=lambda : None,
                sgRNA_info_file=lambda : None,
                uninduced_atc_file=lambda : None,
            )
            def visulize_combined_counts(*args):
                the_file_path = args[-1]
                from pytransit.components.samples_area import samples
                from pytransit.generic_tools import misc
                df = pd.read_csv(the_file_path, sep="\t", comment="#")
                summary_entries = []
                for each in df.columns:
                    if each == "sgRNA":
                        continue
                    samples.wig_table.add({
                        # add hidden link to object
                        "__wig_obj":df[each],
                        # NOTE: all of these names are used by other parts of the code (caution when removing or renaming them)
                        "name": each,
                        "count": len(df),
                        "sum": df[each].values.sum(),
                        "average": round(df[each].values.mean()),
                        "max": round(df[each].values.max()),
                    })
            
            def visulize_metadata(*args):
                metadata_file = args[-1]
                from pytransit.components.samples_area import samples
                from pytransit.generic_tools import misc
                print(f'''args = {args}''')
                df = pd.read_csv(metadata_file, sep="\t", comment="#")
                conditions = misc.no_duplicates(df["column_name"].values.tolist())
                other_columns = [ column for column in df.columns if column != "column_name" ]
                for each in conditions:
                    samples.conditions_table.add({
                        "Condition": each,
                        **{
                            other_column: ", ".join(
                               f"{each}" for each in misc.no_duplicates(df[df["column_name"] == each][other_column].values.tolist()) 
                            )
                                for other_column in other_columns
                        },
                    })
                
                # condition_choices = [ str(each) for each in misc.no_duplicates(df["column_name"].values.tolist())       ] 
                drug_choices =      [ str(each) for each in misc.no_duplicates(df["drug"].values.tolist())              ] 
                days_choices =      [ str(each) for each in misc.no_duplicates(df["days_predepletion"].values.tolist()) ] 
                
                # NOTE: this dynamically adds more textboxes
                self.value_getters.control_condition = panel_helpers.create_choice_input(panel, main_sizer, label="Control condition", options=drug_choices, default_option=None, tooltip_text="")
                self.value_getters.drug = panel_helpers.create_choice_input(panel, main_sizer, label="Drug", options=drug_choices, default_option=None, tooltip_text="")
                self.value_getters.days = panel_helpers.create_choice_input(panel, main_sizer, label="Days", options=days_choices, default_option=None, tooltip_text="")
                panel_helpers.NewPanel.recent.refresh()
            
            def load_folder(*args):
                folder_path = args[-1]
                visulize_combined_counts(f"{folder_path}/combined_counts.txt")
                visulize_metadata(f"{folder_path}/samples_metadata.txt")
            
            if debugging_enabled:
                self.value_getters.cgi_folder            = panel_helpers.create_folder_input(panel, main_sizer, button_label="CGI Folder", tooltip_text="", popup_title="", default_folder=None, default_folder_name="", allowed_extensions='*', after_select=load_folder)
            self.value_getters.combined_counts_file      = panel_helpers.create_file_input(panel, main_sizer, button_label="Combined counts file", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*', after_select=visulize_combined_counts)
            self.value_getters.metadata_file             = panel_helpers.create_file_input(panel, main_sizer, button_label="Metadata file", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*', after_select=visulize_metadata)
            self.value_getters.sgRNA_info_file      = panel_helpers.create_file_input(panel, main_sizer, button_label="sgRNA efficacy file", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*')
            self.value_getters.uninduced_atc_file        = panel_helpers.create_file_input(panel, main_sizer, button_label="uninduced ATC file", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*')
            panel_helpers.create_run_button(
                panel, main_sizer, from_gui_function=self.from_gui
            )
            

    @staticmethod
    def from_gui(frame):
        import pandas as pd
        arguments = LazyDict()
        #
        # call all GUI getters, puts results into respective arguments key-value
        #
        for each_key, each_getter in Method.value_getters.items():
            try:
                arguments[each_key] = each_getter()
            except Exception as error:
                logging.error(
                    f"""Failed to get value of "{each_key}" from GUI:\n{error}"""
                )

        #
        # ask for output path(s)
        #
        arguments.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}_output.tsv",
            output_extensions="Common output extensions (*.tsv,*.dat,*.out)|*.txt;*.tsv;*.dat;*.out;|\nAll files (*.*)|*.*",
        )
        # if user didn't select an output path
        if not arguments.output_path:
            return None
        
        with transit_tools.TimerAndOutputs(
            method_name=f"{Method.identifier} extracting abundance",
            output_paths=[arguments.output_path],
            disable=False,
        ) as timer:
            df = Method.extract_abund(
                combined_counts_file=arguments.combined_counts_file or f"{arguments.cgi_folder}/combined_counts.txt",
                metadata_file=arguments.metadata_file or f"{arguments.cgi_folder}/samples_metadata.txt",
                sgRNA_info_file=arguments.sgRNA_info_file or f"{arguments.cgi_folder}/sgRNA_info.txt",
                uninduced_atc_file=arguments.uninduced_atc_file or f"{arguments.cgi_folder}/uninduced_ATC_counts.txt",
                efficacy_col="",
                orf_col="",
                control_condition=arguments.control_condition,
                drug=arguments.drug,
                days=arguments.days,
                fractional_abundance_file=None,
            )
            frac_abund_file = misc.inject_path_extension(path=arguments.output_path, extension="abundances")
            df.to_csv(frac_abund_file, sep="\t", encoding='utf-8')
        
            print(f'''df = {df}''')
            Method.run_model(
                frac_abund_file=frac_abund_file,
                cdr_output_file=arguments.output_path,
            )

@transit_tools.ResultsFile
class CgiResult:
    column_names = ["Position","Reads","Genes"] 
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
                    title=Method.identifier,
                    heading=self.comments_string or misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[],
                ).Show(),
                "Display Gene": lambda *args: (
                    Method.visualize(
                        fractional_abundances_file=self.values_for_result_table.get("frac_abund_file", None),
                        gene=gui_tools.ask_for_text(label_text="Gene name(s) or ORF Id"),
                        fig_location=None,
                        fixed="",
                        origx=False,
                        origy=False,
                    )
                ),
            })
        )
        
        # 
        # read in data
        # 
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        self.values_for_result_table.update(self.extra_data.get("parameters", {}))
        
        # 
        # get summary stats
        #
        self.values_for_result_table.update({
            # HANDLE_THIS (additional summary_info for results table)
            # examples:
                # f"Gene Count": len(self.rows),
                # f"Adj P Value < {Method.significance_threshold}": len([
                #     1 for each in self.rows
                #         if each.get("Adj P Value", 0) < Method.significance_threshold 
                # ]),
        })
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

