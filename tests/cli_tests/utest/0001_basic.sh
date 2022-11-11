#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/Subu_KO_metadata.txt"
comwig="./src/pytransit/data/Subu_KO_combined_wig.txt"
wig1="src/pytransit/data/glycerol_H37Rv_rep1.wig"
wig2="src/pytransit/data/cholesterol_H37Rv_rep1.wig"


# <combined-wig-path> <annotation .prot_table or GFF3> <metadata path> <condition name for control group> <condition name for experimental group> <output file> [Optional Arguments]
python3 ./src/transit.py "utest" \
    "$comwig" \
    "$annotation" \
    "$metadata" \
    "H37Rv_day0" \
    "H37Rv_day32" \
    "$result_file"
