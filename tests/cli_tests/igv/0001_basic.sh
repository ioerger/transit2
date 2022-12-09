#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/subu_ko.transit/metadata.tsv"
comwig="./src/pytransit/data/subu_ko.transit/comwig.tsv"
wig="src/pytransit/data/cholesterol_glycerol.transit/glycerol_rep1.wig"

python3 ./src/transit.py "export" igv \
    "$wig" \
    "$annotation" \
    "$result_file"
