#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
DATA="src/pytransit/data/"

# can also do GSEA method (on any system of functional categories)
transit pathway_enrichment "tests/cli_tests/resampling/0001_basic.sh.1.result" $DATA/H37Rv_sanger_roles.dat $DATA/sanger_roles.dat "$result_file" -M GSEA