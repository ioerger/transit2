#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/iron.transit/metadata.tsv"
comwig="./src/pytransit/data/iron.transit/comwig.tsv"
DATA="src/pytransit/data/"

# Ontologizer is a specialized method for GO terms
python3 ./src/transit.py pathway_enrichment "tests/cli_tests/resampling/0001_basic.sh.1.result" $DATA/H37Rv_GO_terms.txt $DATA/gene_ontology.1_2.3-11-18.obo "$result_file" --M ONT
