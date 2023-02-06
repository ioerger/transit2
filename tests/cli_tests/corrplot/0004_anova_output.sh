#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/iron.transit/annotation.H37Rv.prot_table"
metadata="./src/pytransit/data/iron.transit/metadata.tsv"
comwig="./src/pytransit/data/iron.transit/comwig.tsv"

python3 ./src/transit.py corrplot "$comwig" "$metadata" "$annotation" ./tests/cli_tests/anova/0001_basic.sh.1.result "$result_file.png" -anova

