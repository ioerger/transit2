#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/Subu_KO_metadata.txt"
comwig="./src/pytransit/data/Subu_KO_combined_wig.txt"
wig="src/pytransit/data/glycerol_H37Rv_rep1.wig"

python3 ./src/transit.py "export" igv \
    "$wig" \
    "$annotation" \
    "$result_file"
