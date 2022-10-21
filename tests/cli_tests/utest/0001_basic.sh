#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/Subu_KO_metadata.txt"
comwig="./src/pytransit/data/Subu_KO_combined_wig.txt"
wig1="src/pytransit/data/glycerol_H37Rv_rep1.wig"
wig2="src/pytransit/data/cholesterol_H37Rv_rep1.wig"

#  <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
python3 ./src/transit.py "utest" \
    "$wig1" \
    "$wig2" \
    "$annotation" \
    "$result_file"
