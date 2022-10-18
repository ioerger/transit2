#!/usr/bin/env bash

#gi <combined_wig> <samples_metadata> <conditionA1> <conditionB1> <conditionA2> <conditionB2> <prot_table> <output_file> 

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/Subu_KO_metadata.txt"
comwig="./src/pytransit/data/Subu_KO_combined_wig.txt"

python3 ./src/transit.py gi \
    "$comwig" \
    "$metadata" \
    "H37Rv_day0"\
    "H37Rv_day32" \
    "Rv2680_day0" \
    "Rv2680_day32" \
    "$annotation" \
    "$result_file"


