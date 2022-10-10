#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"

python3 ./src/transit.py ttnfitness \
    "./src/pytransit/data/glycerol_H37Rv_rep1.wig,./src/pytransit/data/glycerol_H37Rv_rep2.wig" \
    ./src/pytransit/genomes/H37Rv.prot_table \
    ./src/pytransit/genomes/H37Rv.fna \
    ./tests/cli_tests/gumbel/0001_basic.sh.1.result \
    "$result_file.genes" \
    "$result_file.sites"