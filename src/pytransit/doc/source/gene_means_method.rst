.. _gene_means:

Gene Means
===========

It is often useful to reduce a combined_wig file (with insertion counts at each TA site)
to a summary of mean counts for each *gene* in each condition,
which can be done with the 'gene_means' command (in console mode),
or by the following sequence of menu items 
in the GUI: 'Pre-processing'->'Generate'->'Gene Means'.

By default, the gene means are calculated for each individual sample
(column) in the combined_wig file.  However, if desired, the user can
request that counts be averaged over replicates of each condition (as
defined in the metadata file), by including the '-cond' flag. When
averaging replicates together, gene_means uses TTR normalization by
default.



::

 Usage:
    > python3 /pacific/home/ioerger/transit/src/transit.py  gene_means <combined_wig> <metadata_file> <annotation_file> <output_file> [Optional Arguments]

    Optional Arguments:
    --n <string> :=  Normalization method. Default: --n TTR
    --iN <N>     :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: --iN 0
    --iC <N>     :=  Ignore TAs within given percentage (e.g. 5) of C terminus. Default: --iC 0
    -cond        :=  Averages counts over replicates of each condition
    


|

.. rst-class:: transit_sectionend
----
