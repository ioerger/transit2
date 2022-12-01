#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/mc2_155_tamu.prot_table"
metadata="./src/pytransit/data/harman_randomized.metadata.tsv"
comwig="./src/pytransit/data/harman_randomized.comwig.tsv"

annotation="./src/pytransit/data/mc2_155_tamu.prot_table"

# Strains: "DKO" "rel" "WT" "x5849"
# Time Slices: 0days, 7days, 14days, 28days

# asking about "KO_REL" condition
python3.9 ./src/transit.py zinb \
    "$comwig" \
    "$annotation" \
    "$metadata" \
    "$result_file" \
    -include-conditions WT_Starved_d28,rel_Starved_d28,5849_Starved_d28,DKO_Starved_d28 \
    -condition "KO_rel" \
    -interactions KO_5849 \
    -gene uvrA