#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/Subu_KO_metadata.txt"
comwig="./src/pytransit/data/Subu_KO_combined_wig.txt"
wig1="src/pytransit/data/glycerol_H37Rv_rep1.wig"
wig2="src/pytransit/data/cholesterol_H37Rv_rep1.wig"


#  <combined-wig-path> <annotation .prot_table or GFF3> <comma-separated wig_ids for control group> <comma-separated wig_ids experimental group> <output file> [Optional Arguments]
python3 ./src/transit.py "utest" \
    "$comwig" \
    "$metadata" \
    "$annotation" \
    "H37Rv_day32_rep1,H37Rv_day32_rep2" \
    "H37Rv_day0_rep1,H37Rv_day0_rep2" \
    "$result_file"
