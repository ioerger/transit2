#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/iron.transit/metadata.tsv"
comwig="./src/pytransit/data/iron.transit/comwig.tsv"

python3 ./src/transit.py heatmap \
    ./tests/cli_tests/anova/0001_basic.sh.1.result \
    "$result_file.png" \
    -anova
    # --topk <int> \
    # --qval <float> \
    # --low_mean_filter <int> \