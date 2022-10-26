#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result.png"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"

# <combined_wig> <annotation_file> <output.png> [-avg_by_conditions <metadata_file>]
python3 ./src/transit.py corrplot "$comwig" "$annotation" "$result_file" -avg_by_conditions "$metadata"

