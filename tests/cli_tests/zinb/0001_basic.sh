#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/222_samples_metadata_cg.txt"
comwig="./src/pytransit/data/111_cholesterol_glycerol_combined.cwig"

# FIXME: https://stackoverflow.com/questions/64002936/error-in-prettynum-internalformatx-trim-digits-nsmall-width-3l-invalid
python3 ./src/transit.py zinb \
    "$comwig" \
    "$annotation" \
    "$metadata" \
    "$result_file" \
    -group-by "Condition" 
    # -interactions "Cholesterol"
