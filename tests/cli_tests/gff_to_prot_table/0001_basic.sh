#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"
gff_path="./src/pytransit/data/NC_003062.2.gff3"

python3 ./src/transit.py convert gff_to_prot \
    "$gff_path" \
    "$result_file"