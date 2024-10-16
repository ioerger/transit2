# Change log

All notable changes to this project will be documented in this file.


## Version 1.1.5 (2024-10-16)
#### Transit:

Minor changes:
 - ensure that font file is included in PyPi package


## Version 1.1.4 (2024-10-14)
#### Transit:

Minor changes:
 - fixed a bug in TrackView that was caused by deprecated function in recent version of Python Image Library (PIL v10.0)


## Version 1.1.3 (2024-09-14)
#### Transit:

Minor changes:
 - allow wig file pathnames to have spaces in combined_wig and metadata files (e.g. for ANOVA and ZINB)
 - changed default alg for BWA from 'mem' back to 'aln' (see documentation on TPP)


## Version 1.1.2, 2024-04-16
#### TRANSIT2:
  - minor bug fix in 'cgi visualize'

Minor changes:
## Version 1.1.1, 2024-04-12
#### TRANSIT2:

Minor changes:
  - updates to CGI functions (CRISPRi-DR) and documentation
  - 'cgi extract_counts' can now process *.fastq.gz files (automatically unzips them)
  - minor changes to how sgRNA ids are handled
  - minor updates to format of input files like sgRNA_info.txt
  - minor changes to command-line args and flags

## Version 1.1.0, 2024-04-12
#### TRANSIT2:

Minor changes:
  - added empirical Bayes FDR analysis to filter significant interacting genes in CGI
  - added '-no_uninduced' flag to CGI command (see documentation)

## Version 1.0.14, 2024-04-04
#### TRANSIT2:

Minor changes:
  - fixed a bug that caused GUI to hang after popup windows were closed

## Version 1.0.13, 2024-04-04
#### TRANSIT2:

Minor changes:
  - fixed a bug for making combined_wig files in GUI using gff files

## Version 1.0.12, 2024-03-17
#### TRANSIT2:

Minor changes:
  - fix minor bug in CRISPRi-DR (cgi) output

## Version 1.0.11, 2024-03-17
#### TRANSIT2:

Minor changes:
  - minor updates to CRISPRi-DR (cgi) method and documentation


## Version 1.0.10, 2024-02-17
#### TRANSIT2:

Minor changes:
  - fixed .tolist() bug in betageom normalization in 'export combined_wig'


## Version 1.0.9, 2024-01-31
#### TRANSIT2:

Minor changes:
  - a few updates to CRISRPi-DR (CGI)


## Version 1.0.8, 2023-12-20
#### TRANSIT2:

Major changes:
  - added CRISPRi-DR method for analyzing CGI data (Chemical-Genetic Interactions)
    - includes GUI interface
    - see documentation
  - added confidence scores to HMM output

Minor changes:
  - fixed LOESS plots to show genome positional bias before and after correction


## Version 1.0.7, 2023-10-19
#### TRANSIT2:

Minor changes:
  - bug fix in GUI for resampling


## Version 1.0.6 2023-10-18
#### TRANSIT2:

Minor changes:
  - bug fix for corrplot (remove dependence on rpy2)
  - minor edits to documentation


## Version 1.0.5 2023-10-15
#### TRANSIT2:

Minor changes:
  - fix something in .readthedocs.taml


## Version 1.0.4 2023-10-15
#### TRANSIT2:

Minor changes:
  - fix docs on readthedocs by adding config file

	
## Version 1.0.3 2023-10-13
#### TRANSIT2:

Minor changes:
  - Minor bugfix in ttnfitness

	
## Version 1.0.2 2023-07-15
#### TRANSIT2:

Minor fix:
  - ensured that all sub-directories are included in distribution


	
## Version 1.0.1 2023-07-11
#### TRANSIT2:

Updated documentation:
  - clarify that this is Transit2
  - put a link to the original version of Transit
  - update installation instructions for pip and git


## Version 1.0.0 2023-05-31
#### TRANSIT2:

Major new release.
  - Re-implmentation from scratch.
  - Updated GUI and command-line interface.

