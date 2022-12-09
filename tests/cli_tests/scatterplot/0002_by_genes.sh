#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/cholesterol_glycerol.transit/metadata.tsv"
comwig="./src/pytransit/data/cholesterol_glycerol.transit/comwig.tsv"
wig_id1="c1"
wig_id2="c2"

python3 ./src/transit.py scatterplot \
    "$comwig" \
    "$annotation" \
    "$metadata" \
    -samp "$wig_id1,$wig_id2" \
    "$result_file.png" \
    --genes \
    --log