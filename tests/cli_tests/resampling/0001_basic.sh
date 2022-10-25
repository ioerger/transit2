#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv_dev.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"

python3 ./src/transit.py resampling \
    -s 1000 \
    -c ./src/pytransit/data/cholesterol_glycerol_combined.dat \
    ./src/pytransit/data/samples_metadata_cg.txt \
    Glycerol \
    Cholesterol \
    ./src/pytransit/data/genomes/H37Rv.prot_table \
    "$result_file"