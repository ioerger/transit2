#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="src/pytransit/data/cholesterol_glycerol.transit/metadata.tsv"
comwig="src/pytransit/data/cholesterol_glycerol.transit/comwig.tsv"

# FIXME: https://stackoverflow.com/questions/64002936/error-in-prettynum-internalformatx-trim-digits-nsmall-width-3l-invalid
python3 ./src/transit.py zinb \
    "$comwig" \
    "$annotation" \
    "$metadata" \
    "$result_file" \
    -group-by "Condition" 
    # -interactions "Cholesterol"
