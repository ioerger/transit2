#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"
condition1="Cholesterol"
condition2="Glycerol"

python3 ./src/transit.py scatterplot \
    "$comwig" \
    "$annotation" \
    "$metadata" \
    -cond \
    "$condition1,$condition2" \
    "$result_file.png" \
    --genes 