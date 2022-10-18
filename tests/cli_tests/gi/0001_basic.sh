#!/usr/bin/env bash

#use Rv2680 Knockout data to get 4 datasets
#gi <combined_wig> <samples_metadata> <conditionA1> <conditionB1> <conditionA2> <conditionB2> <prot_table> <output_file> 

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/samples_metadata_cg.txt"
comwig="./src/pytransit/data/cholesterol_glycerol_combined.dat"

python3 ./src/transit.py gi \
    "$comwig" \
    "$metadata" \
    #conditionA1
    #conditionA2
    #conditionB1
    #conditionB2
    "$annotation" \
    "$result_file"


