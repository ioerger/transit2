#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
DATA="src/pytransit/data/"

# Ontologizer is a specialized method for GO terms
transit pathway_enrichment "tests/cli_tests/resampling/0001_basic.sh.1.result" $DATA/H37Rv_GO_terms.txt $DATA/GO_term_names.dat "$result_file" -M ONT
