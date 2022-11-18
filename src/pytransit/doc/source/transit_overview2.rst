



TRANSIT Overview (4.0)
======================


The analysis methods in Transit are also described in this `PDF manual
<https://orca1.tamu.edu/essentiality/transit/transit-manual.pdf>`_ .

overview

TPP
---

* explain how to run TPP

* outputs wig files (coord of TA sites, ins cnt); ref genome at top

* look at .tnseq_stats files for density, etc.

* Tn5 - prefix, etc; all coords in wigs


Combined_wig and Metadata Files
-------------------------------

This is a new step; multiple wig files are combined into a single CW file.
metadata describes each sample.  Most commonly, several replicates will be 
associated with a specific condition.
The rest of the analyses in Transit have been
redesigned to take CW and meta files as inputs, which simplifies things.

* export combined_wig command (must be same genome)

* required columns in metadata (make with Excel, save as tab-sep)


Genome Annotations (prot_tables and gff files)
------------------

* explain what a prot_table is and why we prefer them

* for analysis, it is critical to use the right prot_table corresponding to the ref genome
used to make wigs with TPP (why?)

* can convert from gff to prot_table

* prot_table adjustment?


Command Line
------------

flags



GUI
---

* loading files

* sample actions

  *  sample dropdowns tasks ("select tool" - loess, track-view...)

* param panel for methods

* display table, analysis-specific actions 

* message bar (errors)


Pre-Processing
--------------

* tnseq_stats, QQplots, track-view, scatterplots, corrplot, gene_means

* QC - data quality, diagnostics

* normalization (TTR, betageom)

* Analyses (for Himar1 datasets)

 * 3 types:

 * single

 * pairwise 

 * multiple

 * output files (tab-sep spreadsheets)

 * hits are ususally Qval<0.05

* Analysis for Tn5

Results and Post-Processing
---------------------------

* display table

* volcano plot

* Pathway Enrichment Analysis





Developers
----------

=======================  ============  ==============================================================================
 Name                    Time Active          Contact Information
=======================  ============  ==============================================================================
Thomas R. Ioerger        2015-Present  `http://faculty.cs.tamu.edu/ioerger/ <http://faculty.cs.tamu.edu/ioerger/>`_
Michael A. DeJesus       2015-2018     `http://mad-lab.org <http://mad-lab.org>`_
Chaitra Ambadipudi       2015
Eric Nelson              2016
Siddharth Subramaniyam   2018
Sanjee Choudhery         2021-
Jeff Hykin               2022-
=======================  ============  ==============================================================================




References
----------


If you use TRANSIT, please cite the following reference:


.. [DeJesus2015TRANSIT] `DeJesus, M.A., Ambadipudi, C., Baker, R., Sassetti, C., and Ioerger, T.R. (2015). TRANSIT - a Software Tool for Himar1 TnSeq Analysis. PLOS Computational Biology, 11(10):e1004401 <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004401>`_



Development of TRANSIT is funded by the National Institutes of Health (www.nih.gov/) grant U19 AI107774.



Other references for methods utilized by TRANSIT:



.. [DeJesus2013]  `DeJesus, M.A., Zhang, Y.J., Sassettti, C.M., Rubin, E.J.,
  Sacchettini, J.C., and Ioerger, T.R. (2013). Bayesian analysis of gene essentiality based on sequencing of transposon insertion libraries. Bioinformatics, 29(6):695-703. <http://www.ncbi.nlm.nih.gov/pubmed/23361328>`_


.. [DeJesus2013HMM] `DeJesus, M.A., Ioerger, T.R. A Hidden Markov Model for identifying essential and growth-defect regions in bacterial genomes from transposon insertion sequencing data. BMC Bioinformatics. 2013. 14:303 <http://www.ncbi.nlm.nih.gov/pubmed/24103077>`_


.. [DeJesus2014] `DeJesus, M.A. and Ioerger, T.R. (2014). Capturing uncertainty by modeling local transposon insertion frequencies improves discrimination of essential genes. IEEE Transactions on Computational Biology and Bioinformatics, 12(1):92-102. <http://www.ncbi.nlm.nih.gov/pubmed/26357081>`_



.. [DeJesus2016] `DeJesus, M.A. and Ioerger, T.R. (2016). Normalization of transposon-mutant library sequencing datasets to improve identification of conditionally essential genes. Journal of Bioinformatics and Computational Biology, 14(3):1642004 <http://www.ncbi.nlm.nih.gov/pubmed/26932272>`_


.. [DeJesus2017NAR] `DeJesus, M.A., Nambi, S., Smith, C.M., Baker, R.E., Sassetti, C.M., Ioerger, T.R. Statistical analysis of genetic interactions in Tn-Seq data.  Nucleic Acids Research. 2017. 45(11):e93. doi: 10.1093/nar/gkx128. <https://www.ncbi.nlm.nih.gov/pubmed/28334803>`_

.. [ZINB] `Subramaniyam S, DeJesus MA, Zaveri A, Smith CM, Baker RE, Ehrt S, Schnappinger D, Sassetti CM, Ioerger TR. (2019).  Statistical analysis of variability in TnSeq data across conditions using Zero-Inflated Negative Binomial regression. *BMC Bioinformatics*. 2019 Nov 21;20(1):603. doi: 10.1186/s12859-019-3156-z. <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3156-z>`_

.. [Choudhery2021] `Choudhery S, Brown AJ, Akusobi C, Rubin EJ, Sassetti CM, Ioerger TR. Modeling Site-Specific Nucleotide Biases Affecting Himar1 Transposon Insertion Frequencies in TnSeq Data Sets. *mSystems*. 2021 Oct 26;6(5):e0087621. doi: 10.1128/mSystems.00876-21. <https://pubmed.ncbi.nlm.nih.gov/34665010/>`_
