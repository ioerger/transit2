
.. _normalization:

Normalization
=============


Proper normalization is important as it ensures that other sources of
variability are not mistakenly treated as real differences in
datasets. TRANSIT provides various normalization methods, which are
briefly described below:

- **TTR:**
    Trimmed Total Reads (TTR), normalized by the total
    read-counts (like totreads), but trims top and bottom 5% of
    read-counts. **This is the recommended normalization method for most cases**
    as it has the beneffit of normalizing for difference in
    saturation in the context of resampling.

- **nzmean:**
    Normalizes datasets to have the same mean over the
    non-zero sites.

- **totreads:**
    Normalizes datasets by total read-counts, and scales
    them to have the same mean over all counts.

- **zinfnb:**
    Fits a zero-inflated negative binomial model, and then
    divides read-counts by the mean. The zero-inflated negative
    binomial model will treat some empty sites as belonging to the
    "true" negative binomial distribution responsible for read-counts
    while treating the others as "essential" (and thus not influencing
    its parameters).

- **quantile:**
    Normalizes datasets using the quantile normalization
    method described by `Bolstad et al.
    (2003) <http://www.ncbi.nlm.nih.gov/pubmed/12538238>`_. In this
    normalization procedure, datasets are sorted, an empirical
    distribution is estimated as the mean across the sorted datasets
    at each site, and then the original (unsorted) datasets are
    assigned values from the empirical distribution based on their
    quantiles.

- **betageom:**
    Normalizes the datasets to fit an "ideal" Geometric
    distribution with a variable probability parameter *p*. Specially
    useful for datasets that contain a large skew. See :ref:`BGC` .

- **nonorm:**
    No normalization is performed.

Command-line
------------

You can call Transit to normalize wig files or combined_wig files from the command-line.
(It will automatically determine the type of input file.)

::

 Usage:
     > python3 transit.py normalize <wig_file or combined_wig_file> <output_file> [Optional Arguments]

     Optional Arguments:
     --n <string>  :=  Normalization method (e.g. TTR, betageom, etc). Default: --n TTR



GUI
---

Normalization may be applied to combined_wig files loaded into the GUI 
by selecting the following sequence of menu items:
'Pre-processing'->'Generate'->'A normalized combined_wig file using...'


Examples
-------

::

  > python3 src/transit.py normalize Rv_1_H37RvRef.wig Rv_1_H37RvRef_TTR.wig -n TTR

  > python3 src/transit.py normalize Rv_1_H37RvRef.wig Rv_1_H37RvRef_BG.wig -n betageom



.. rst-class:: transit_sectionend
----
