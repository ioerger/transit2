#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/mc2_155_tamu.prot_table"
metadata="./src/pytransit/data/randomized_harman.transit/metadata.tsv"
comwig="./src/pytransit/data/randomized_harman.transit/comwig.tsv"

# Strains: "DKO" "rel" "WT" "x5849"
# Time Slices: 0days, 7days, 14days, 28days

# asking about "KO_REL" condition
python3 ./src/transit.py zinb \
    "$comwig" \
    "$metadata" \
    "$annotation" \
    "$result_file" \
    --include-conditions WT_Starved_d28,rel_Starved_d28,5849_Starved_d28,DKO_Starved_d28 \
    --condition "KO_rel" \
    --interactions KO_5849 \
    --gene uvrA