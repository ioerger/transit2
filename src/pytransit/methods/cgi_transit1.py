import os
import time
import math
import random
import numpy
import scipy.stats
import datetime
import pandas

from pytransit.analysis import base
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools

from rpy2.robjects import (
    r,
    DataFrame,
    globalenv,
    IntVector,
    FloatVector,
    StrVector,
    packages as rpackages,
)

############# Description ##################

short_name = "CGI"
long_name = "Chemical Genetic Analysis"
short_desc = "CGI Analysis of CRISPRi libraries"
long_desc = "CGI Analysis of CRISPRi libraries"
transposons = []

columns = ["Position", "Reads", "Genes"]

class CGI_Method(base.SingleConditionMethod):
    def __init__(self):
        ctrldata = None  # initializers for superclass
        annotation_path = ""
        output_file = ""
        replicates = "Sum"
        normalization = "nonorm"
        LOESS = False
        ignoreCodon = True
        NTerminus = 0.0
        CTerminus = 0.0
        wxobj = None
        # this initialization seems pointless for CGI, but must do this for base class...
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
            NTerminus=NTerminus,
            CTerminus=CTerminus,
            wxobj=wxobj,
        )

    @classmethod
    def usage_string(self):
        return """usage (6 sub-commands):
    python3 ../src/transit.py CGI extract_counts <fastq file> <ids file> <counts file>
    python3 ../src/transit.py CGI create_combined_counts <comma seperated headers> <counts file 1> <counts file 2> ... <counts file n>  <combined counts file>
    python3 ../src/transit.py CGI extract_abund <combined counts file> <metadata file> <control condition> <sgRNA strength file> <uninduced ATC file> <drug> <days>  <fractional abundundance file>
    python3 ../src/transit.py CGI run_model <fractional abundundance file>  <CRISPRi DR results file>
    python3 ../src/transit.py CGI visualize <fractional abundance> <gene> <output figure location>
        -fixed xmin=x,xmax=x,ymin=y,ymax=y : set the values you would to be fixed in this comma seperated format. Not all values need to be set for ex, a valid arguement is "xmin=0,ymax=5"
        -origx : flag to turn on original scale axes rather than log scale for Concentration default=off
        -origy : flag to turn on original scale axes rather than log scale for Realtive Abundances default=off
"""

    @classmethod
    def fromargs(self, rawargs):
        if not has_r:
            print("Error: R and rpy2 (~= 3.0) required to run heatmap.")
            print(
                "After installing R, you can install rpy2 using the command \"pip install 'rpy2~=3.0'\""
            )
            sys.exit(0)

        (args, kwargs) = transit_tools.cleanargs(rawargs)
        if len(args) < 1:
            print(self.usage_string())
        self.cmd = args[0]
        self.args = args[1:]
        self.kwargs = kwargs

        return self()

    def Run(self):
        cmd, args, kwargs = self.cmd, self.args, self.kwargs

        if cmd == "extract_counts":
            if len(args) < 2:
                print("You have provided incorrect number of args")
                print(self.usage_string())
                sys.exit(0)
            fastq_file = args[0]
            ids_file = args[1]
            counts_file = args[2]
            self.extract_counts(fastq_file, ids_file, counts_file)

        elif cmd == "create_combined_counts":
            if len(args) < 2:
                print("You have provided incorrect number of args")
                print(self.usage_string())
                sys.exit(0)
            headers = args[0].split(",")
            counts_file_list = args[1:-1]
            combined_counts = args[-1]
            self.create_combined_counts(headers, counts_file_list, combined_counts)

        elif cmd == "extract_abund":
            if len(args) < 7:
                print("You have provided incorrect number of args")
                print(self.usage_string())
                sys.exit(0)
            combined_counts_file = args[0]
            metadata_file = args[1]
            control_condition = args[2]
            sgRNA_strength_file = args[3]
            no_dep_abund = args[4]
            drug = args[5]
            days = args[6]
            fractional_abundance_file = args[7]
            self.extract_abund(
                combined_counts_file,
                metadata_file,
                control_condition,
                sgRNA_strength_file,
                no_dep_abund,
                drug,
                days,
                fractional_abundance_file,
            )
        elif cmd == "run_model":
            if len(args) < 1:
                print("You have provided incorrect number of args")
                print(self.usage_string())
                sys.exit(0)
            ifile_path = args[0]  # example frac_abund_RIF_D5.txt
            cdr_output_file = args[1]
            self.run_model(ifile_path, cdr_output_file)
        elif cmd == "visualize":
            if len(args) < 3:
                print("You have provided incorrect number of args")
                print(self.usage_string())
                sys.exit(0)
            frac_abund_file = args[0]
            gene = args[1]
            fig_location = args[2]

            flags = "-fixed -origx -origy".split()
            for arg in args:
                if arg[0] == "-" and arg not in flags:
                    print("flag unrecognized: %s" % arg)
                    print(self.usage_string())
                    sys.exit(0)

            fixed_axes = kwargs.get(
                "fixed", ""
            )  # did it this way rather than having xmin, xmax etc. as seperate flag because the - in say -2.5 would trigger flag recognition
            origx_axis = kwargs.get("origx", False)
            origy_axis = kwargs.get("origy", False)

            self.visualize(
                frac_abund_file,
                gene,
                fig_location,
                fixed=fixed_axes,
                origx=origx_axis,
                origy=origy_axis,
            )

        else:
            print("You have not entered a valid command, here are the options")
            print(self.usage_string())

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

        # A,B = "AGCTTCTTTCGAGTACAAAAAC","TCCCAGATTATATCTATCACTGA"
        A, B = "GTACAAAAAC", "TCCCAGATTA"
        lenA = len(A)
        cnt, nreads, recognized = 0, 0, 0
        for line in open(fastq_file):
            cnt += 1
            if cnt % 4 == 2:
                nreads += 1
                if nreads % 1000000 == 0:
                    sys.stderr.write(
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
            sys.stderr.write("Adding in file # %s \n" % f)
            counts_df = pd.read_csv(f, sep="\t")
            counts_df["sgRNA"] = counts_df[counts_df.columns[0]].str.split(
                "_v", expand=True
            )[0]
            counts_df = counts_df.drop(columns=[counts_df.columns[0]])
            counts_df.set_index("sgRNA", inplace=True)
            df_list.append(counts_df)
        combined_df = pd.concat(df_list, axis=1)
        combined_df.columns = headers
        combined_df.to_csv(combined_counts_file, sep="\t")
        sys.stderr.write(
            "Number of sgRNAs in combined counts file (present in all counts files): %d \n"
            % len(combined_df)
        )

    def extract_abund(
        self,
        combined_counts_file,
        metadata_file,
        control_condition,
        sgRNA_strength_file,
        no_dep_abund,
        drug,
        days,
        fractional_abundance_file,
        PC=1e-8,
    ):
        import pandas as pd

        metadata = pd.read_csv(metadata_file, sep="\t")
        metadata = metadata[
            ((metadata["drug"] == drug) | (metadata["drug"] == control_condition))
            & (metadata["days_predepletion"] == int(days))
        ]
        if len(metadata) == 0:
            sys.stderr.write(
                "This combination of conditions does not exist in your metadata file. Please select one that does"
            )
            sys.exit(0)
        elif drug not in metadata["drug"].values.tolist():
            sys.stderr.write(
                "%s is not found in your metadata. Add the drug's information in the metadata file or select a different drug"
                % drug
            )
            sys.exit(0)
        elif int(days) not in metadata["days_predepletion"].values.tolist():
            sys.stderr.write(
                "%d is not found in your metadata days of predepletion column. Add the day's information in the metadata file or select a different day"
                % days
            )
            sys.exit(0)
        elif control_condition not in metadata["drug"].values.tolist():
            sys.stderr.write(
                "%s is not found in your metadata. Add the corresponding information in the metadata file or select a different control"
                % control_condition
            )
            sys.exit(0)
        metadata = metadata.sort_values(by=["conc_xMIC"])
        column_names = metadata["column_name"].values.tolist()
        concs_list = metadata["conc_xMIC"].values.tolist()

        print("# Condition Tested : " + drug + " D" + days)
        headers = []
        combined_counts_df = pd.read_csv(combined_counts_file, sep="\t", index_col=0)
        combined_counts_df = combined_counts_df[column_names]

        if len(combined_counts_df.columns) == 0:
            sys.stderr.write(
                "The samples assocaited with the selected drugs do not exist in your combined counts file. Please select one that does and check your metadata file has corresponding column names"
            )
            sys.exit(0)
        elif len(combined_counts_df.columns) < len(metadata):
            sys.stderr.write(
                "WARNING: Not all of the samples from the metadata based on this criteron have a column in the combined counts file"
            )

        sgRNA_strength = pd.read_csv(sgRNA_strength_file, sep="\t", index_col=0)
        sgRNA_strength = sgRNA_strength.iloc[:, -1:]
        sgRNA_strength.columns = ["sgRNA strength"]
        sgRNA_strength["sgRNA"] = sgRNA_strength.index
        sgRNA_strength["sgRNA"] = sgRNA_strength["sgRNA"].str.split("_v", expand=True)[
            0
        ]
        sgRNA_strength.set_index("sgRNA", inplace=True)

        no_dep_df = pd.read_csv(no_dep_abund, sep="\t", index_col=0, header=None)
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
        sys.stderr.write("Disregarding Empty or Negative sgRNAs\n")
        sys.stderr.write(
            "%d sgRNAs are all of the following files : sgRNA strength metadata, uninduced ATC counts file, combined counts file\n"
            % len(abund_df)
        )
        f = open(fractional_abundance_file, "w")
        headers = ["sgRNA strength", "uninduced ATC values"]
        for i, col in enumerate(column_names):
            abund_df[col] = abund_df[col] / abund_df[col].sum()
            abund_df[col] = (abund_df[col] + PC) / (
                abund_df["uninduced ATC values"] + PC
            )
            headers.append(str(concs_list[i]) + "_" + str(i))
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

        abund_df.to_csv(f, sep="\t", index=False)

    #####################################################

    # derived from logsigmoidfit.R
    # see heatmap.py for example of how to put data in a pandas.DataFrame and call an R function like make_heatmapFunc()

    def run_model(self, frac_abund_file, cdr_output_file):
        import pandas as pd
        import numpy as np
        from mne.stats import fdr_correction
        import statsmodels.api as sm

        frac_abund_df = pd.read_csv(frac_abund_file, sep="\t", comment="#")

        drug_output = []
        for i, gene in enumerate(set(frac_abund_df["gene"])):
            # print(i,gene)
            sys.stderr.write("Analyzing Gene # %d \n" % i)
            gene_df = frac_abund_df[frac_abund_df["gene"] == gene]
            orf = gene_df["orf"].iloc[0]
            gene_df = gene_df.drop(columns=["orf", "gene", "uninduced ATC values"])

            melted_df = gene_df.melt(
                id_vars=["sgRNA", "sgRNA strength"], var_name="conc", value_name="abund"
            )
            melted_df["conc"] = (
                melted_df["conc"].str.split("_", expand=True)[0].astype(float)
            )
            min_conc = min(melted_df[melted_df["conc"] > 0]["conc"])
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
            sys.stderr.flush()

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
        sys.stderr.write("%d Total Significant Gene Interactions\n" % n)
        sys.stderr.write("%d Significant Gene Depletions\n" % depl_n)
        sys.stderr.write("%d Significant Gene Enrichments\n" % enrich_n)

        drug_out_df = drug_out_df.replace(r"\s+", np.nan, regex=True).replace(
            "", np.nan
        )
        drug_out_df.to_csv(cdr_output_file, sep="\t", index=False)

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

        abund_df = pd.read_csv(fractional_abundances_file, sep="\t", comment="#")
        with open(fractional_abundances_file) as f:
            first_line = f.readline()
            condition = first_line.split(" : ")[1]

        abund_df = abund_df[(abund_df["gene"] == gene) | (abund_df["orf"] == gene)]
        if len(abund_df) == 0:
            sys.stderr.write("Gene not found : %d \n" % idx)
            sys.exit(0)
        abund_df = abund_df.reset_index(drop=True)
        all_slopes = []

        df_list = []
        for idx, row in abund_df.iterrows():
            sys.stderr.write("Fitting sgRNA # : %d \n" % idx)
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

        plot_df = pd.concat(df_list)

        plt.figure()
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
            print(
                "You have provided more than one xmin value. Figure will be created using flexible xmin"
            )
        elif len(xmin_list) == 1:
            xmin_temp = float(xmin_list[0].split("=")[1])
            if xmin_temp > xmax:
                print(
                    "You have provided an xmin value greater than the xmax value. Figure will be created using flexible xmin"
                )
            else:
                xmin = xmin_temp

        xmax_list = [x for x in range_values if re.search("xmax=", x)]
        if len(xmax_list) > 1:
            print(
                "You have provided more than one xmax value. Figure will be created using flexible xmax"
            )
        elif len(xmax_list) == 1:
            xmax_temp = float(xmax_list[0].split("=")[1])
            if xmax_temp < xmin:
                print(
                    "You have provided an xmax value less than the xmin value. Figure will be created using flexible xmax"
                )
            else:
                xmax = xmax_temp

        ymin_list = [x for x in range_values if re.search("ymin=", x)]
        if len(ymin_list) > 1:
            print(
                "You have provided more than one ymin value. Figure will be created using flexible ymin"
            )
        elif len(ymin_list) == 1:
            ymin_temp = float(ymin_list[0].split("=")[1])
            if ymin_temp < ymax:
                print(
                    "You have provided an ymin value greater than the ymax value. Figure will be created using flexible ymin"
                )
            else:
                ymin = ymin_temp

        ymax_list = [x for x in range_values if re.search("ymax=", x)]
        if len(ymax_list) > 1:
            print(
                "You have provided more than one ymax value. Figure will be created using flexible ymax"
            )
        elif len(ymax_list) == 1:
            ymax_temp = float(ymax_list[0].split("=")[1])
            if ymax_temp < ymin:
                print(
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
        plt.gca().set_title(gene + "\n" + condition, wrap=True)
        plt.tight_layout()
        plt.savefig(fig_location)


################################

if __name__ == "__main__":
    G = CGI.fromargs(sys.argv[1:])
    G.Run()
