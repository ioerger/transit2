.. _tnseq_stats:

.. rst-class:: 
Tnseq_stats
===========

You can generate the same table to statistics as on the Quality Control panel in the GUI
from the command-line using the 'tnseq_stats' command.  Here is an example:

::

  > python3 src/transit.py tnseq_stats -help

  usage: python3 src/transit.py tnseq_stats <file.wig>+ [-o <output_file>]
         python3 src/transit.py tnseq_stats -c <combined_wig> [-o <output_file>]

  > python3 src/transit.py tnseq_stats -c src/pytransit/data/cholesterol_glycerol.transit/comwig.tsv

  Dataset                    Density Mean_ct NZmean NZmedian Max_ct   Total_cts Skewness Kurtosis Pickands_Tail_Index
  cholesterol_H37Rv_rep1.wig 0.439   139.6   317.6  147      125355.5 10414005   54.8     4237.7  0.973
  cholesterol_H37Rv_rep2.wig 0.439   171.4   390.5  148      704662.8 12786637  105.8    14216.2  1.529
  cholesterol_H37Rv_rep3.wig 0.359   173.8   484.2  171      292294.8 12968502   42.2     2328.0  1.584
  glycerol_H37Rv_rep1.wig    0.419   123.3   294.5  160        8813.3  9195672    4.0      33.0   0.184
  glycerol_H37Rv_rep2.wig    0.516   123.8   240.1  127        8542.5  9235984    4.0      33.5   0.152


The output file is tab-separated text file (spreadsheet) with the following columns:

+----------------------+-----------------------------------------------------------------+
| Column Header        | Column Definition                                               |
+======================+=================================================================+
| Dataset              | name of sample (.wig file)                                      |
+----------------------+-----------------------------------------------------------------+
| Density              | saturation (percent of TA sites with non-zero insertion counts  |
+----------------------+-----------------------------------------------------------------+
| Mean_ct              | the mean count over all TA sites                                |
+----------------------+-----------------------------------------------------------------+
| NZmean               | the mean count over non-zero TA sites                           |
+----------------------+-----------------------------------------------------------------+
| NZmedian             | the median count over non-zero TA sites                         |
+----------------------+-----------------------------------------------------------------+
| Max_ct               | highest count over all TA sites (to check for outliers)         |
+----------------------+-----------------------------------------------------------------+
| Total_cts            | total insertion counts summed over all TA sites                 |
+----------------------+-----------------------------------------------------------------+
| Skewness             | 3rd-order moment of read count distribution                     |
+----------------------+-----------------------------------------------------------------+
| Kurtosis             | 4th-order moment of read count distribution                     |
+----------------------+-----------------------------------------------------------------+
| Pickands Tail Index  | another measure of skewness of read count distributions         |
+----------------------+-----------------------------------------------------------------+


**Signs of potential problems with a dataset:**

 * low density (<30%)
 * low NZmean (<10) - note: this can be affected by normalization of input combined_wig file (might need to re-run 'export combined_wig' with '\-\-n nonorm' to see means of raw data)
 * max_ct: in most datasets, this is usually in the range of thousands to tens of thousands; if it is over a million, this could be an outlier (super-high counts at one or a few sites, possibly due to positive selection), which could be throwing the rest of the counts off
 * skewness: it is difficult to give a hard cutoff, but skewness > 50 could be a sign that a sample is noisy
 * Pickands Tail Index: it is difficult to give a hard cutoff, but PTI>0.5 could be a sign that a sample is noisy (skewed), and PTI>1 is bad.

Pickands Tail Index (PTI) is defined in: `James Pickands III. (1975) 
"Statistical Inference Using Extreme Order Statistics."
Ann. Statist. 3(1):119-131. <https://doi.org/10.1214/aos/1176343003>`_
It is calculated using a formula based on the order statistics of the
distribution of counts (highest counts in sorted order), and increases for
distributions with heavier tails (outliers).  In Transit, the PTI is
calculated over the highest counts, with ranks 10-100.

The analysis methods in Transit will work with noisy samples, but the
results could be affected (e.g. reduced sensitivity of detecting
conditionally essential genes).  Two options to consider are: 1)
dropping the noisiest samples, 2) applying a non-linear normalization
like the Beta-Geometric Correction (BGC) ('\-\-n betagenom').

See also :ref:`Quality Control <transit_quality_control>` for a discussion of
how to interpret these metrics, and what to do if you have noisy samples.



.. **To Do: (TRI, 12/19/22)**

..  * explain PTI
..  * check for PTI in output files and usage()
..  * check headers in output file (does it match what is here?)
..  * add PTI to transit-quality-control (in transit_features.rst)


.. rst-class:: transit_sectionend
----
