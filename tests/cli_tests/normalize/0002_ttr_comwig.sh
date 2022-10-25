#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"
wig="./src/pytransit/data/111_cholesterol_glycerol_combined.cwig"

python3 ./src/transit.py normalize \
    -c "$comwig"\
    "$result_file" \
    -n TTR