#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"

# check if file exists (macrophages isnt in the git repo)
if [ -f "src/pytransit/data.ignore/macrophages.comwig" ]
then
    python3 ./src/transit.py anova \
        ./src/pytransit/data.ignore/macrophages.comwig \
        ./src/pytransit/data.ignore/macrophages_metadata.txt \
        ./src/pytransit/genomes/H37Rv.prot_table \
        "$result_file" \
        --ref Untreated \
        --exclude-conditions Input \
        -n nonorm
fi