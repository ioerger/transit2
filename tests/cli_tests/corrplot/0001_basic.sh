#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"

transit corrplot ./src/pytransit/data/cholesterol_glycerol_combined.dat \ "$result_file"