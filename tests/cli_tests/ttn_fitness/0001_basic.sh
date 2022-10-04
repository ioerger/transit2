#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"

python3 ./src/transit.py ttnfitness \
    "./src/pytransit/data/glycerol_H37Rv_rep1.wig,./src/pytransit/data/glycerol_H37Rv_rep2.wig" \
    ./src/pytransit/genomes/H37Rv.prot_table \
    ./src/pytransit/genomes/H37Rv.fna \
    ./tests/data/gumbel.1.output \
    "$result_file.genes" \
    "$result_file.sites"