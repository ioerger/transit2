#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"

python3 ./src/transit.py hmm \
    "src/pytransit/data/cholesterol_H37Rv_rep3.wig,src/pytransit/data/cholesterol_H37Rv_rep2.wig,src/pytransit/data/cholesterol_H37Rv_rep1.wig" \
    ./src/pytransit/data/genomes/H37Rv.prot_table \
    "$result_file"