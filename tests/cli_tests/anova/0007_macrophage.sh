#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"

# check if file exists (macrophages isnt in the git repo)
if [ -f "src/pytransit/data.ignore/macrophages.comwig" ]
then
    python3 ./src/transit.py anova \
        ./src/pytransit/data.ignore/macrophages.comwig \
        ./src/pytransit/data.ignore/macrophages_metadata.txt \
        ./src/pytransit/data/genomes/H37Rv.prot_table \
        "$result_file" \
        -ref Untreated \
        -exclude-conditions Input \
        -n nonorm
fi