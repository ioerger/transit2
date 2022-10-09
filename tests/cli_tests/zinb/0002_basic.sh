#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"

python3 ./src/transit.py zinb \
    "./src/pytransit/data/111_cholesterol_glycerol_combined.cwig" \
    "./src/pytransit/data/222_samples_metadata_cg.txt" \
    "$result_file" \
    --condition "Glycerol" \
    --interactions "Cholesterol"