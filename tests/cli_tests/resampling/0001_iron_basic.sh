#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/iron.transit/metadata.tsv"
comwig="./src/pytransit/data/iron.transit/comwig.tsv"
condition1="HighFeMBT"
condition2="LowFeMBT"

# <combined_wig_file> <annotation_file> <metadata_file> <ctrl_condition> <exp_condition> <output_file> [Optional Arguments]
python3 ./src/transit.py resampling \
    --s 1000 \
    ./src/pytransit/data/iron.transit/comwig.tsv \
    ./src/pytransit/data/iron.transit/metadata.tsv \
    ./src/pytransit/data/genomes/H37Rv.prot_table \
    "$condition1" \
    "$condition2" \
    "$result_file"