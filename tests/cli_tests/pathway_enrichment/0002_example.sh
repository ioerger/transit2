#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/cholesterol_glycerol.transit/metadata.tsv"
comwig="./src/pytransit/data/cholesterol_glycerol.transit/comwig.tsv"
DATA="src/pytransit/data/"

# with new COG 20 categories
python3 ./src/transit.py pathway_enrichment "tests/cli_tests/resampling/0001_basic.sh.1.result" $DATA/H37Rv_COG_20_roles.txt $DATA/COG_20_roles.txt "$result_file"