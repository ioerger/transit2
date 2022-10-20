#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"

if python3 ./src/transit.py anova \
    ./src/pytransit/data/cholesterol_glycerol_combined.dat \
    ./src/pytransit/data/samples_metadata_cg.txt \
    ./src/pytransit/data/genomes/H37Rv.prot_table \
    "$result_file" \
    -alpha 100 \
    --ref Untreated
then
    false
else
    true # we expect/want this case to fail (checking helpful error message)
fi