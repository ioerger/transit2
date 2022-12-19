#!/usr/bin/env bash


#gi <combined_wig> <samples_metadata> <conditionA1> <conditionB1> <conditionA2> <conditionB2> <prot_table> <output_file> 

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/subu_ko.transit/metadata.tsv"
comwig="./src/pytransit/data/subu_ko.transit/comwig.tsv"

python3 ./src/transit.py gi \
    "$comwig" \
    "$annotation" \
    "$metadata" \
    "H37Rv_day0"\
    "H37Rv_day32" \
    "Rv2680_day0" \
    "Rv2680_day32" \
    "$result_file"


