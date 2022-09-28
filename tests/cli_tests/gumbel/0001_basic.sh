#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"

python3 ./src/transit.py gumbel \
    "/Users/jeffhykin/repos/transit/src/pytransit/data/glycerol_H37Rv_rep1.wig,/Users/jeffhykin/repos/transit/src/pytransit/data/glycerol_H37Rv_rep2.wig" \
    ./src/pytransit/genomes/H37Rv.prot_table \
    "$result_file"