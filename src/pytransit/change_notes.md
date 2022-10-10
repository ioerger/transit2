naming check:
    - 'TA Sites', 'Coordinates', 'Sites'
    - "Gene Name"
    - "LFC", "Log FC", "Log 2 FC"
    - "Non Insertions"
    - "GD Sites", "GD Counts", "GD" in HMM
    - "State Call", "Essentiality Call", "Call"
    - "ES Count", "ES GD NE GA N/A", "ES Sites", "ES"


TTN Fitness
    SitesFile:
        comlumn_names = [
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
    GenesFile:
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
    strings that are possibly keys
        "A"
        "Actual LFC"
        "Analysis"
        "b"
        "bit"
        "C"
        "c"
        "Condition"
        "Coordinates"
        "Count Of Sites With Insertions"
        "Description"
        "Display Table"
        "Downstream TTN"
        "E"
        "EB"
        "ES"
        "ESB"
        "Essentiality Call"
        "Filtering ES/ESB Genes"
        "Filtering Short Genes"
        "Fitness Ratio"
        "G"
        "g"
        "GA"
        "GD"
        "Gene Plus TTN M1 Adj P Value"
        "Gene Plus TTN M1 Coef"
        "Gene Plus TTN States"
        "Gene Saturation"
        "Genes"
        "Getting Data"
        "Getting Genome"
        "glycerol_gumbel.out"
        "himar1"
        "Id"
        "igr"
        "index"
        "Insertion Count"
        "k"
        "Local Average"
        "m"
        "M1 Adj P Value"
        "M1 Coef"
        "M1 P Value"
        "M1 Pred Log Count"
        "M1 Predicted Count"
        "Making Fitness Estimations"
        "Mean Insertion Count"
        "Gene Name"
        "NE"
        "nonorm"
        "Normalizing using: %s"
        "Number Of Insertions Within ORF"
        "ORF"
        "processing data"
        "r"
        "Resampling - Volcano plot"
        "Sites"
        "State"
        "T"
        "TA"
        "Total TA Site Count"
        "TTN Fitness Assessment"
        "TTN Fitness"
        "TTNFitness Summary"
        "ttnfitness.test"
        "TTR"
        "U"
        "Upstream TTN"
        "Writing File: %s"
        "Writing To Output Files"
        "Writing To Output"
        "y"

Gumbel
    column_names = [
        "ORF", 
        "Gene Name", 
        "Description", 
        "Number Of Insertions Within ORF",
        "Total Number Of TA Sites Within ORF", 
        "Length Of Maximum Run Of Non Insertions",
        "Nucleotide Span For Maximum Run Of Non Insertions",
        "Posterior Probability Of Essentiality", # Z Bar
        "Essentiality Call",
    ]
    strings that are possibly keys
        "b"
        "Burnin"
        "Condition to analyze:"
        "Condition"
        "Description"
        "Display Table"
        "Doing Regression"
        "E"
        "EB"
        "Essentiality Call"
        "Getting Data from %s"
        "Getting Data"
        "Gumbel"
        "himar1"
        "iC"
        "Id"
        "If the density of the dataset is too low, the Gumbel method will not work."
        "iN"
        "Length Of Maximum Run Of Non Insertions"
        "n"
        "Gene Name"
        "NE"
        "nonorm"
        "Normalizing using: %s"
        "Nucleotide Span For Maximum Run Of Non Insertions"
        "Number Of Insertions Within ORF"
        "ORF"
        "Posterior Probability Of Essentiality"
        "r"
        "s"
        "S"
        "Samples"
        "Sum"
        "t"
        "This is likely to have been caused by poor data (e.g. too sparse)."
        "Total Number Of TA Sites Within ORF"
        "trim"
        "TTR"
        "U"

HMM 
    column_names = [
        "ORF",
        "Gene",
        "Annotation",
        "TAs",
        "ES Sites",
        "GD Sites",
        "NE Sites",
        "GA Sites",
        "Saturation",
        "Mean",
        "Call",
    ]
    strings that are possibly keys
        "Annotation"
        "Calculating pins"
        "Call"
        "cli"
        "Combining Replicates as '%s'"
        "Correct for Genome Positional Bias"
        "Creating HMM Genes Level Output"
        "Description"
        "Display Table"
        "ES Count"
        "ES GD NE GA N/A"
        "ES Sites"
        "ES"
        "Essentiality Call"
        "GA Count"
        "GA Sites"
        "GA"
        "GD Count"
        "GD Sites"
        "GD"
        "Gene Count"
        "Gene"
        "gui"
        "Hidden Markov Model"
        "himar1"
        "HMM_Genes"
        "HMM_Sites"
        "HMM"
        "iC"
        "ignore"
        "iN"
        "l"
        "Length Of Maximum Run Of Non Insertions"
        "Location"
        "Mean Insertions"
        "Mean Reads"
        "Mean"
        "n"
        "N/A"
        "Gene Name"
        "NE Count"
        "NE Sites"
        "NE"
        "nonorm"
        "Normalizing using: %s"
        "Nucleotide Span For Maximum Run Of Non Insertions"
        "Number Of Insertions Within ORF"
        "ORF"
        "Padj"
        "Padj<{Analysis.significance_threshold}"
        "parameters"
        "Performing loess_correction Correction"
        "Please select a condition"
        "Posterior Probability Of Essentiality"
        "Probability ES"
        "Probability GA"
        "Probability GD"
        "Probability NE"
        "r"
        "Read Count"
        "Replicates:"
        "Saturation"
        "State Call"
        "State"
        "Sum"
        "TAs"
        "Total Number Of TA Sites Within ORF"
        "Total Sites"
        "TTR"
        "Unknown State"
        "w"

Pathway Enrichment strings that are possibly keys
    'phypergeometric'
    "Adj P Value"
    "Analysis"
    "COG_20"
    "Description"
    "Display Table"
    "Enrichment Exponent"
    "Enrichment"
    "Expected"
    "FET"
    "Gene Count"
    "Genes In Path"
    "Genes"
    "GSEA"
    "H37Rv-COG"
    "H37Rv-GO"
    "H37Rv-Sanger"
    "himar1"
    "K Plus PC"
    "LFC_col"
    "LFC"
    "Log 2 FC"
    "log2FC"
    "M"
    "Method"
    "name:"
    "Not a valid method"
    "Nperm"
    "Number Adjusted By PC"
    "Number of Permutations"
    "ONT"
    "P Value"
    "p-value"
    "p"
    "Padj"
    "Padj<{Analysis.significance_threshold}"
    "parameters"
    "Pathway Enrichment"
    "Pathway"
    "PC"
    "Pval_col"
    "Qval_col"
    "Ranking"
    "ranking"
    "Significant Genes"
    "Significent Genes In Path"
    "SLPV"
    "Smeg-COG"
    "Smeg-GO"
    "SPLV"
    "tn5"
    "Total Genes"

Anova Columns Names & Accessed Column Names
    "Rv",
    "Gene Name",
    "TAs",
    *[ f"Mean {condition_name}" for condition_name in conditions_list ],
    *[  f"LFC {condition_name}" for condition_name in conditions_list ],
    "MSR",
    "MSE With Alpha",
    "Fstat",
    "P Value",
    "Adj P Value",
    "Status"
    
Resampling Column Names & Accessed Columns Names
    "ORF",
    "Gene Name",
    "Description",
    "Sites",
    "Mean Control",
    "Mean Experimental",
    "Log 2 FC",
    "Sum Control",
    "Sum Experimental",
    "Delta Mean",
    "P Value",
    "Z Score", # sometimes excluded
    "Adj P Value",

GI Column Names & Accessed Columns Names 
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

Tnseq Stats Column Names & Accessed Columns Names 
    "Dataset",
    "Density",
    "Mean Count",
    "Non Zero Mean",
    "Non Zero Median",
    "Max Count",
    "Total Counts",
    "Skewness",
    "Kurtosis",
    "Pickands Tail Index",

Combined Wig Columns
    "TA Site Position"
    "ORF",
    "Gene Name"


some replacements that have already happened as CLI commands (for Jeff)
```shell
sd '"Orf"' '"ORF"' src/**/*.py
sd 'from pytransit.universal_data *import universal'   'from pytransit.interfaces import gui, cli'                              src/**/*.py
sd 'universal\.interface == "gui"'                     'gui.is_active'                                                          src/**/*.py
sd 'universal\.interface != "gui"'                     'not gui.is_active'                                                      src/**/*.py
sd 'universal\.interface != "cli"'                     'gui.is_active'                                                          src/**/*.py
sd 'universal\.interface == "cli"'                     'not gui.is_active'                                                      src/**/*.py
sd 'universal\.'                                       'gui.'                                                                   src/**/*.py
sd 'from pytransit.interfaces *import gui, cli'        'from pytransit.globals import gui, cli, root_folder, debugging_enabled' src/**/*.py
sd 'gui.debugging_enabled'                             'debugging_enabled'                                                      src/**/*.py
sd 'gui.root_folder'                                   'root_folder'                                                            src/**/*.py
sd 'gui\.interface *== *("|'"'"')gui("|'"'"')'         'gui.is_active'                                                          src/**/*.py
sd 'gui\.interface *!= *("|'"'"')gui("|'"'"')'         'not gui.is_active'                                                      src/**/*.py
sd 'gui\.interface *!= *("|'"'"')console("|'"'"')'     'gui.is_active'                                                          src/**/*.py
sd 'gui\.interface *== *("|'"'"')console("|'"'"')'     'not gui.is_active'                                                      src/**/*.py
sd 'FIX ME'                                            'FIXME'                                                                  src/**/*.py
sd 'sys.argv'                                          'console_tools.subcommand_prefix#FIXME'                                  src/methods/**/*.py
sd '\bAnalysis\b'                                      'Method'                                                                 src/methods/**/*.py
sd 'Method = GUI = Method'                             ''                                                                       src/**/*.py
```