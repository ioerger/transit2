#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/mc2_155_tamu.prot_table"
metadata="./src/pytransit/data/Harman_metadata_v3.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"


# Strains: "DKO" "rel" "WT" "x5849"
# Time Slices: 0days, 7days, 14days, 28days


# What is a valid input for --interactions
# What is a valid input for --covars
# Are they all boolean? no
# Are numbers treated differently? yes

# asking about "KO_REL" condition
python3 ./src/transit.py zinb \
    "$combined_wig_metadata_v3" \
    "$metadata" \
    "$annotation" \
    "$result_file" \
    --include-conditions WT_Starved_d28,rel_Starved_d28,5849_Starved_d28,DKO_Starved_d28 \
    --group-by KO_rel \
    --interactions KO_5849

# > python3 ../../transit/src/transit.py zinb Harman_combined_wig_smeg_BGC.txt combined_wig_metadata_v3.txt mc2_155_tamu.prot_table ZINB_GI_DKO_BGC_d28.txt --include-conditions WT_Starved_d0,rel_Starved_d0,5849_Starved_d0,DKO_Starved_d0,WT_Starved_d28,rel_Starved_d28,5849_Starved_d28,DKO_Starved_d28 --condition Strain --interactions Starvation --prot_table mc2_155_tamu.prot_table


# python3 ./src/transit.py zinb \
#     "./src/pytransit/data/111_cholesterol_glycerol_combined.cwig" \
#     "./src/pytransit/data/222_samples_metadata_cg.txt" \
#     "$result_file" \
#     --interactions "Glycerol"