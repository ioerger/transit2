#!/usr/bin/env bash

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"
annotation="./src/pytransit/data/genomes/H37Rv.prot_table"
metadata="./src/pytransit/data/iron.transit/metadata.tsv"
comwig="./src/pytransit/data/iron.transit/comwig.tsv"
wig_id1="HighFeMBT"
wig_id2="HighFeMBT2"
wig_id3="HighFeMBT3"
wig_id4="LowFeMBT"
wig_id5="LowFeMBT2"
wig_id6="FeCMBT"
wig_id7="FeCMBT2"
wig_id8="Hemin"
wig_id9="Hemin2"
wig_id10="Hemin3"
wig_id11="Hemoglobin"
wig_id12="Hemoglobin2"
wig_id13="HeminMBT"
wig_id14="HeminMBT2"

python3 ./src/transit.py scatterplot \
    "$comwig" \
    "$metadata" \
    "$annotation" \
    --samp "$wig_id1,$wig_id2" \
    "$result_file.png"