#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="src/pytransit/data/iron.transit/metadata.tsv"
comwig="src/pytransit/data/iron.transit/comwig.tsv"

python3 ./src/transit.py normalize \
    "$comwig" \
    "$result_file" \
    --n TTR