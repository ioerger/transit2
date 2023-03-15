#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/iron.transit/annotation.H37Rv.prot_table"
# annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/iron.transit/metadata.tsv"
comwig="./src/pytransit/data/iron.transit/comwig.tsv"
condition1="HighFeMBT"
condition2="LowFeMBT"
# HighFeMBT
# LowFeMBT
# FeCMBT
# Hemin
# Hemoglobin
# HeminMBT


# <combined_wig_file> <annotation_file> <metadata_file> <ctrl_condition> <exp_condition> <output_file> [Optional Arguments]
python3 ./src/transit.py resampling \
    --s 1000 \
    -a \
    "$comwig" \
    "$metadata" \
    "$annotation" \
    "$condition1" \
    "$condition2" \
    "$result_file"