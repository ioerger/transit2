#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/cholesterol_glycerol.transit/metadata.tsv"
comwig="./src/pytransit/data/cholesterol_glycerol.transit/comwig.tsv"

python3 ./src/transit.py hmm \
    "src/pytransit/data/cholesterol_glycerol.transit/cholesterol_rep3.wig,src/pytransit/data/cholesterol_glycerol.transit/cholesterol_rep2.wig,src/pytransit/data/cholesterol_glycerol.transit/cholesterol_rep1.wig" \
    ./src/pytransit/data/genomes/H37Rv.prot_table \
    "$result_file"