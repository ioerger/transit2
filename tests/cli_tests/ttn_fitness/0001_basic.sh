#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"

python3 ./src/transit.py ttnfitness_gui \
    "/Users/jeffhykin/repos/transit/src/pytransit/data/glycerol_H37Rv_rep1.wig,/Users/jeffhykin/repos/transit/src/pytransit/data/glycerol_H37Rv_rep2.wig" \
    ./src/pytransit/genomes/H37Rv.prot_table \
    H37Rv.fna \
    "$result_file" \
    "$result_file.genes.out" \
    "$result_file.sites.out"