#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
DATA="src/pytransit/data/"

# with COG categories
transit pathway_enrichment "tests/cli_tests/resampling/0001_basic.sh.1.result" $DATA/H37Rv_COG_roles.dat $DATA/COG_roles.dat "$result_file"