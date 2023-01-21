#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/cholesterol_glycerol.transit/annotation.H37Rv.prot_table"
# annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/cholesterol_glycerol.transit/metadata.tsv"
comwig="./src/pytransit/data/cholesterol_glycerol.transit/comwig.tsv"
wig1="src/pytransit/data/cholesterol_glycerol.transit/glycerol_rep1.wig"
wig2="src/pytransit/data/cholesterol_glycerol.transit/cholesterol_rep1.wig"


# <combined_wig_file> <annotation_file> <metadata_file> <ctrl_condition> <exp_condition> <output_file> [Optional Arguments]
python3 ./src/transit.py resampling \
    -s 1000 \
    --a \
    "$comwig" \
    "$metadata" \
    "$annotation" \
    Glycerol \
    Cholesterol \
    "$result_file"