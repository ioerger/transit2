
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
(a unusually high insertion count at just a single TA site).
A reference for this method is `(Santa Maria et al., 2014)
<https://www.ncbi.nlm.nih.gov/pubmed/25104751>`__.

Usage
-----


::

  > python3 /pacific/home/ioerger/transit/src/transit.py  utest <combined-wig-path> <annotation .prot_table or GFF3> <metadata path> <condition name for control group> <condition name for experimental group> <output file> [Optional Arguments]

  Optional Arguments:
  -n <string>     :=  Normalization method. Default: -n TTR
  -iz             :=  Include TA sites with zero across conditions.
  -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Turned Off.
  -iN <float>     :=  Ignore TAs occuring at given fraction (as integer) of the N terminus. Default: -iN 0
  -iC <float>     :=  Ignore TAs occuring at given fraction (as integer) of the C terminus. Default: -iC 0


[Note to Jeff: we should probably get rid of -l, since LOESS is not relevant for U-test.]

[Note to Jeff: isn't -iz true by default?]

[Note to Jeff: normally I would includes a section explaining the parameters, but they are so obvious in this case.]

[Note to Jeff: update usage for combined wigs]


Output
------

Explain the columns in the output file...

Point out that the significant genes are those with Padj<0.05.

Runtime
-------

If this typically takes less than a minute, say so...


|

.. rst-class:: transit_sectionend
----
