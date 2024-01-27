#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result.png"
wig="./tests/data/test.wig"

python3 ./src/transit.py loess \
    "$wig" \
    "$result_file"