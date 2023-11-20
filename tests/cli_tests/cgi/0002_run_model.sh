#!/usr/bin/env bash


#gi <combined_wig> <samples_metadata> <conditionA1> <conditionB1> <conditionA2> <conditionB2> <prot_table> <output_file> 

result_file="./tests/cli_tests/$(basename "$(dirname "$0")")/$(basename "$0").1.result"

ids_file="data/CGI/IDs.H37Rv.CRISPRi.lib.txt"
combined_counts_file="data/CGI/RIF_D1_combined_counts.txt"
metadata_file="data/CGI/samples_metadata.txt"
uninduced_atc_file="data/CGI/uninduced_ATC_counts.txt"
sg_rna_strength_file="data/CGI/sgRNA_info.txt"
fractional_abundance_file="data/CGI/fractional_abundance_file.txt"

drug="RIF"
days="1"
control_condition="RIF"
# fastq_file
# ids_file
# headers_comma_separated
# output_counts_file
# uninduced_atc_file
# gene
# crispri_dr_results_file
 
# python3 ./src/transit.py cgi extract_counts "$fastq_file" "$ids_file" "$output_counts_file"
# python3 ./src/transit.py cgi create_combined_counts "$headers_comma_separated" "$counts_file1" "$counts_file1" "$combined_counts_file"
# python3 ./src/transit.py cgi extract_abund "$combined_counts_file" "$metadata_file" "$control_condition" "$sg_rna_strength_file" "$uninduced_atc_file" "$drug" "$days" "$result_file"
python3 ./src/transit.py cgi run_model "$fractional_abundance_file" "$result_file"
# python3 ./src/transit.py cgi visualize "$fractional_abundance_file" "$gene" "$output_figure_location"

# Optional Arguments: 
#     -fixed xmin=x,xmax=x,ymin=y,ymax=y := set the values you would to be fixed in this comma seperated format. Not all values need to be set for ex, a valid arguement is "xmin=0,ymax=5"
#     -origx := flag to turn on original scale axes rather than log scale for Concentration default=off
#     -origy := flag to turn on original scale axes rather than log scale for Realtive Abundances default=off
    