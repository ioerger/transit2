import sys
import os
import time
import ntpath
import math
import random
import datetime
import collections
import heapq

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
        "--fixed",  # fixed axes
        "-origx", # non-log axes
        "-origy", # non-log axes
    ]
    usage_string = f"""
        Usage (6 sub-commands):

        {console_tools.subcommand_prefix} {cli_name} extract_counts <fastq file> <ids file> <output counts file>
        
        {console_tools.subcommand_prefix} {cli_name} create_combined_counts <comma seperated headers> <counts file 1> <counts file 2> ... <counts file n> <combined counts file>
        
        {console_tools.subcommand_prefix} {cli_name} extract_abund <combined counts file> <metadata file> <control condition> <sgRNA strength file> <uninduced ATC file> <drug> <days>  <fractional abundundance file>
        
        {console_tools.subcommand_prefix} {cli_name} run_model <fractional abundundance file>  <CRISPRi DR results file>
        
        {console_tools.subcommand_prefix} {cli_name} visualize <fractional abundance> <gene> <output figure location> [Optional Arguments]
        Optional Arguments: 
            -fixed xmin=x,xmax=x,ymin=y,ymax=y := set the values you would to be fixed in this comma seperated format. Not all values need to be set for ex, a valid arguement is "xmin=0,ymax=5"
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
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=3)

        fastq_file= args[1]
        ids_file = args[2]
        counts_file = args[3]
        Method.extract_counts(fastq_file, ids_file, counts_file)

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
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=8)

        combined_counts_file      = args[0]
        metadata_file             = args[1]
        control_condition         = args[2]
        sg_rna_strength_file      = args[3]
        uninduced_atc_file        = args[4]
        drug                      = args[5]
        days                      = args[6]
        fractional_abundance_file = args[7]
        Method.extract_abund(
            combined_counts_file,
            metadata_file,
            control_condition,
            sg_rna_strength_file,
            uninduced_atc_file,
            drug,
            days,
            fractional_abundance_file,
        )

    @staticmethod
    @cli.add_command(cli_name, "run_model")
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=2)

        ifile_path = args[0]
        cdr_output_file = args[1]
        Method.run_model(ifile_path, cdr_output_file)

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

    def extract_counts(self, fastq_file, ids_file, counts_file):
        f = open(counts_file, "w")
        IDs = []
        barcodemap = {}  # hash from barcode to full ids
        for line in open(ids_file):
            w = line.rstrip().split("\t")
            id = w[0]
            v = id.split("_")
            if len(v) < 3:
                continue
            barcode = v[2]
            IDs.append(id)
            # reverse-complement of barcodes appears in reads, so hash them that way
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
        f.close()

    def create_combined_counts(self, headers, counts_list, combined_counts_file):
        import pandas as pd

        df_list = []
        for f in counts_list:
            logging.log("Adding in file # %s \n" % f)
            counts_df = pd.read_csv(f, sep="\t", comment="#")
            counts_df["sgRNA"] = counts_df[counts_df.columns[0]].str.split(
                "_v", expand=True
            )[0]
            counts_df = counts_df.drop(columns=[counts_df.columns[0]])
            counts_df.set_index("sgRNA", inplace=True)
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
        sg_rna_strength_file,
        uninduced_atc_file,
        drug,
        days,
        fractional_abundance_file,
        PC=1e-8,
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

        logging.log("# Condition Tested : " + str(drug) + " D" + str(days))
        headers = []
        combined_counts_df = pd.read_csv(combined_counts_file, sep="\t", index_col=0, comment="#")
        combined_counts_df = combined_counts_df[column_names]

        if len(combined_counts_df.columns) == 0:
            logging.error(
                "The samples assocaited with the selected drugs do not exist in your combined counts file. Please select one that does and check your metadata file has corresponding column names"
            )
        elif len(combined_counts_df.columns) < len(metadata):
            logging.log(
                "WARNING: Not all of the samples from the metadata based on this criteron have a column in the combined counts file"
            )

        sgRNA_strength = pd.read_csv(sg_rna_strength_file, sep="\t", index_col=0, comment="#")
        sgRNA_strength = sgRNA_strength.iloc[:, -1:]
        sgRNA_strength.columns = ["sgRNA strength"]
        sgRNA_strength["sgRNA"] = sgRNA_strength.index
        sgRNA_strength["sgRNA"] = sgRNA_strength["sgRNA"].str.split("_v", expand=True)[
            0
        ]
        sgRNA_strength.set_index("sgRNA", inplace=True)

        no_dep_df = pd.read_csv(uninduced_atc_file, sep="\t", index_col=0, header=None, comment="#")
        no_dep_df = no_dep_df.iloc[:, -1:]
        no_dep_df.columns = ["uninduced ATC values"]
        no_dep_df["uninduced ATC values"] = (
            no_dep_df["uninduced ATC values"] / no_dep_df["uninduced ATC values"].sum()
        )
        no_dep_df["sgRNA"] = no_dep_df.index
        no_dep_df["sgRNA"] = no_dep_df["sgRNA"].str.split("_v", expand=True)[0]
        no_dep_df.set_index("sgRNA", inplace=True)

        abund_df = pd.concat([sgRNA_strength, no_dep_df, combined_counts_df], axis=1)
        abund_df = abund_df[
            ~(
                abund_df.index.str.contains("Negative")
                | abund_df.index.str.contains("Empty")
            )
        ]
        logging.log("Disregarding Empty or Negative sgRNAs")
        logging.log(
            "%d sgRNAs are all of the following files : sgRNA strength metadata, uninduced ATC counts file, combined counts file"
            % len(abund_df)
        )

        headers = ["sgRNA strength", "uninduced ATC values"]
        if fractional_abundance_file != None:
            f = open(fractional_abundance_file, "w")
        for i, col in enumerate(column_names):
            abund_df[col] = abund_df[col] / abund_df[col].sum()
            abund_df[col] = (abund_df[col] + PC) / (
                abund_df["uninduced ATC values"] + PC
            )
            headers.append(str(concs_list[i]) + "_" + str(i))
            if fractional_abundance_file != None:
                f.write("# " + str(concs_list[i]) + " conc_xMIC" + " - " + col + "\n")

        abund_df.columns = headers
        abund_df["sgRNA"] = abund_df.index.values.tolist()
        abund_df[["orf-gene", "remaining"]] = abund_df["sgRNA"].str.split(
            "_", n=1, expand=True
        )
        abund_df[["orf", "gene"]] = abund_df["orf-gene"].str.split(":", expand=True)
        abund_df = abund_df.drop(columns=["orf-gene", "remaining"])
        abund_df = abund_df.dropna()

        abund_df.insert(0, "sgRNA strength", abund_df.pop("sgRNA strength"))
        abund_df.insert(0, "uninduced ATC values", abund_df.pop("uninduced ATC values"))
        abund_df.insert(0, "gene", abund_df.pop("gene"))
        abund_df.insert(0, "orf", abund_df.pop("orf"))
        abund_df.insert(0, "sgRNA", abund_df.pop("sgRNA"))
        
        if fractional_abundance_file != None:
            abund_df.to_csv(f, sep="\t", index=False)
        
        return abund_df

    def run_model(self, frac_abund_file, cdr_output_file, disable_logging=False, frac_abund_df=None):
        from pytransit.components import panel_helpers, parameter_panel
        import pandas as pd
        import numpy as np
        from mne.stats import fdr_correction
        import statsmodels.api as sm
        
        if type(frac_abund_df) == type(None):
            frac_abund_df = pd.read_csv(frac_abund_file, sep="\t", comment="#")
        
        items = set(frac_abund_df["gene"])
        total = len(items)
        drug_output = []
        for progress, gene in ProgressBar(items, title="Analyzing Genes"):
            gene_index = progress.index
            percentage = (100.0 * gene_index / (total-1))
            gene_df = frac_abund_df[frac_abund_df["gene"] == gene]
            orf = gene_df["orf"].iloc[0]
            gene_df = gene_df.drop(columns=["orf", "gene", "uninduced ATC values"])

            melted_df = gene_df.melt(
                id_vars=["sgRNA", "sgRNA strength"], var_name="conc", value_name="abund"
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
                drug_output.append([orf, gene, len(gene_df)] + [np.nan] * 6)
                continue

            Y = melted_df["logsig abund"]
            X = melted_df.drop(columns=["abund", "logsig abund", "sgRNA", "conc"])
            X = sm.add_constant(X)
            model = sm.OLS(Y, X)
            results = model.fit()
            coeffs = results.params
            pvals = results.pvalues
            drug_output.append(
                [orf, gene, len(gene_df)]
                + coeffs.values.tolist()
                + pvals.values.tolist()
            )
            if gui.is_active:
                parameter_panel.progress_update(f"Analyzing Genes  {percentage:.2f}%\r", percentage)

        drug_out_df = pd.DataFrame(
            drug_output,
            columns=[
                "Orf",
                "Gene",
                "Nobs",
                "intercept",
                "coefficient sgRNA_strength",
                "coefficient concentration dependence",
                "pval intercept",
                "pval pred_logFC",
                "pval concentration dependence",
            ],
        )

        mask = np.isfinite(drug_out_df["pval concentration dependence"])
        pval_corrected = np.full(
            drug_out_df["pval concentration dependence"].shape, np.nan
        )
        pval_corrected[mask] = fdr_correction(
            drug_out_df["pval concentration dependence"][mask]
        )[1]
        drug_out_df["qval concentration dependence"] = pval_corrected
        drug_out_df = drug_out_df.replace(np.nan, 1)

        drug_out_df["Z"] = (
            drug_out_df["coefficient concentration dependence"]
            - drug_out_df["coefficient concentration dependence"].mean()
        ) / drug_out_df["coefficient concentration dependence"].std()
        drug_out_df["Significant Interactions"] = [0] * len(drug_out_df)
        drug_out_df.loc[
            (drug_out_df["qval concentration dependence"] < 0.05)
            & (drug_out_df["Z"] < -2),
            "Significant Interactions",
        ] = -1
        drug_out_df.loc[
            (drug_out_df["qval concentration dependence"] < 0.05)
            & (drug_out_df["Z"] > 2),
            "Significant Interactions",
        ] = 1
        drug_out_df.insert(
            0, "Significant Interactions", drug_out_df.pop("Significant Interactions")
        )

        n = len(drug_out_df[drug_out_df["Significant Interactions"] != 0])
        depl_n = len(drug_out_df[drug_out_df["Significant Interactions"] == -1])
        enrich_n = len(drug_out_df[drug_out_df["Significant Interactions"] == 1])
        logging.log("%d Total Significant Gene Interactions" % n)
        logging.log("%d Significant Gene Depletions" % depl_n)
        logging.log("%d Significant Gene Enrichments" % enrich_n)

        drug_out_df = drug_out_df.replace(r"\s+", np.nan, regex=True).replace(
            "", np.nan
        )
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
                ),
            )
        except Exception as error:
            print("error", error)
            import code; code.interact(local={**globals(),**locals()})
        # drug_out_df.to_csv(cdr_output_file, sep="\t", index=False)
    
    # logging.log("Done")

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
        abund_df_lowercased["gene"] = tuple(each.lower() for each in abund_df["gene"].values)
        abund_df_lowercased["orf"] = tuple(each.lower() for each in abund_df["orf"].values)
        gene_no_case = gene.lower()
        abund_df_filtered = abund_df_lowercased[(abund_df_lowercased["gene"] == gene_no_case) | (abund_df_lowercased["orf"] == gene_no_case)]
        if len(abund_df_filtered) == 0:
            available_genes = list(set(abund_df["gene"].tolist()+abund_df["orf"].tolist()))
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
                data["sgRNA strength"] = [row["sgRNA strength"]] * len(data)
                data["slope"] = [results.params[1]] * len(data)
                data["Concentration"] = raw_X
                data["Relative Abundance"] = raw_Y
                df_list.append(data)
            except Exception as error:
                print(error)
                import code; code.interact(local={**globals(),**locals()})

        plot_df = pd.concat(df_list)

        cmap = mpl.colors.LinearSegmentedColormap.from_list(
            "", ["#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500"], N=len(abund_df)
        )
        palette = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]

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
                hue="sgRNA strength",
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
                hue="sgRNA strength",
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
                hue="sgRNA strength",
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
                hue="sgRNA strength",
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
                vmin=plot_df["sgRNA strength"].min(), vmax=0, clip=False
            ),
            cmap=cmap,
        )
        g.figure.colorbar(sm1, shrink=0.8, aspect=50, label="sgRNA strength")
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
                sg_rna_strength_file=lambda : None,
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
            self.value_getters.sg_rna_strength_file      = panel_helpers.create_file_input(panel, main_sizer, button_label="sgRNA strength file", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*')
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
                sg_rna_strength_file=arguments.sg_rna_strength_file or f"{arguments.cgi_folder}/sgRNA_info.txt",
                uninduced_atc_file=arguments.uninduced_atc_file or f"{arguments.cgi_folder}/uninduced_ATC_counts.txt",
                control_condition=arguments.control_condition,
                drug=arguments.drug,
                days=arguments.days,
                fractional_abundance_file=None,
            )
            frac_abund_file = misc.inject_path_extension(path=arguments.output_path, extension="abundances")
            df.to_csv(frac_abund_file, sep="\t", encoding='utf-8', index=False)
        
            print(f'''df = {df}''')
            Method.run_model(
                frac_abund_file=frac_abund_file,
                frac_abund_df=df,
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

