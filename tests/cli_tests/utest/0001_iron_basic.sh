#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/iron.transit/annotation.H37Rv.prot_table"
# annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/iron.transit/metadata.tsv"
comwig="./src/pytransit/data/iron.transit/comwig.tsv"
condition1="HighFeMBT"
condition2="LowFeMBT"

# <combined-wig-path> <annotation .prot_table or GFF3> <metadata path> <condition name for control group> <condition name for experimental group> <output file> [Optional Arguments]
python3 ./src/transit.py "utest" \
    "$comwig" \
    "$metadata" \
    "$annotation" \
    "$condition1" \
    "$condition2" \
    "$result_file"
