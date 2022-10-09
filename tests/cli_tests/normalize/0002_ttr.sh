#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"

python3 ./src/transit.py normalize \
    "./src/pytransit/data/111_cholesterol_glycerol_combined.cwig" \
    ./src/pytransit/genomes/H37Rv.prot_table \
    "$result_file" \
    -n TTR