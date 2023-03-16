.. _Utest:

Mann-Whitney U-test (utest)
===========================

This is a method for comparing datasets
from a TnSeq library evaluated in
two different conditions (i.e. pairwise-comparison, analogous to resampling).
This is a non-parametric *rank-based* test on whether the level of insertions in a
gene or chromosomal region are significantly higher or lower in one
condition than the other.  Effectively, the insertion counts at the TA
sites in the region are pooled and sorted.  Then the combined ranks of the counts
in region A are compared to those in region B, and p-value is calculated
that reflects whether there is a significant difference in the ranks.

The results (genes with adjusted P-value<0.05) should be similar to what
you get with resampling.
The advantage of this method is that it is less sensitive to outliers
(e.g an unusually high insertion count at just a single TA site).
A reference for this method is `(Santa Maria et al., 2014)
<https://www.ncbi.nlm.nih.gov/pubmed/25104751>`__.

Usage
-----


::

  > python3 transit.py utest <combined_wig_file> <metadata_file> <annotation_file> <condition_for_control> <condition_for_experimental> <output_file> [Optional Arguments]

  Optional Arguments:
    --n <string>     :=  Normalization method. Default: --n TTR
    --iN <float>     :=  Ignore TAs occurring at given fraction (as integer) of the N terminus. Default: --iN 0
    --iC <float>     :=  Ignore TAs occurring at given fraction (as integer) of the C terminus. Default: --iC 0


.. [Note to Jeff: we should probably get rid of -l, since LOESS is not relevant for U-test.]

.. [Note to Jeff: isn't -iz true by default?]

.. [Note to Jeff: normally I would includes a section explaining the parameters, but they are so obvious in this case.]

.. [Note to Jeff: change the order of CL args: <combined_wig> <metadata> <annotation> - it should be like this for all methods.  Check the usage() strings and from_args().  Also, can you update the usage blocks in each method in the documentation?]

GUI Mode
-------

Utest can be access though the "Method" tab in the Menu Bar.
    .. image:: _images/utest_selection_gui.png
       :width: 1000
       :align: center 

The parameters to input through the parameter panel for the method is equivalent to the command line usage, except
in the GUI format we name the output files using the prefix passed in.
    .. image:: _images/utest_parameter_panel.png
       :width: 1000
       :align: center

Output
------

The output file is tab-separated text file (spreadsheet) with the following columns:

+-----------------+-----------------------------------------------------------------+
| Column Header   | Column Definition                                               |
+=================+=================================================================+
| Orf             | Gene ID.                                                        |
+-----------------+-----------------------------------------------------------------+
| Name            | Name of the gene.                                               |
+-----------------+-----------------------------------------------------------------+
| Desc            | Gene description.                                               |
+-----------------+-----------------------------------------------------------------+
| Sites           | Number of TA sites in the gene.                                 |
+-----------------+-----------------------------------------------------------------+
| Mean Ctrl       | Mean of read counts in condition 1. (avg over TA sites and reps)|
+-----------------+-----------------------------------------------------------------+
| Mean Exp        | Mean of read counts in condition 2.                             |
+-----------------+-----------------------------------------------------------------+
| Log 2 FC        | Log-fold-change of exp (treatment) over ctrl (untreated)        |
+-----------------+-----------------------------------------------------------------+
| U Statistic     | test statistic reflecting which condition has higher counts     |
+-----------------+-----------------------------------------------------------------+
| P Value         | 2-tailed P-value of u_stat based on Mann-Whitney                |
+-----------------+-----------------------------------------------------------------+
| Adj P Value     | Adjusted p-value controlling for the FDR (Benjamini-Hochberg)   |
+-----------------+-----------------------------------------------------------------+

The U-test statistic is a value reflecting whether insertion counts at TA sites in control samples are higher on average than counts in the experimental samples.

*Important:* The significant genes (conditionally essential) are those with Padj<0.05.

The significant genes should be comparable with those from resampling.


Runtime
-------

The utest method is relatively fast, and should take less than a minute on most datasets.


|

.. rst-class:: transit_sectionend
----
