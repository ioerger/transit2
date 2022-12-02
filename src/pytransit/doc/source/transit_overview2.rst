



TRANSIT Overview (4.0)
======================


Transit is a python-based software system that combines statistical
algorithms for analyzing TnSeq data (from sequencing transposon-mutant
libraries).  It has GUI (graphical-user interface) to make it easy for
users to process their data, conduct essentiality analyses, and visualize results.
Transit routines can also be run from the command line.  It works in
Linux, MacOS, and Windows.  Transit requires various python
packages to be installed, as well as R (a statistical program that Transit calls for
some analyses).

Transit has two main phases: 

 * a pre-processing phase (called TPP) which analyzes raw sequencing data (fastq files) and extracts transposon insertion counts at genomic coordinates (i.e. TA dinucleotides, for the Himar1 transposon), and   
 * statistical analyses (using Transit proper).   

The statistical analyses are focused on identifying genes
that are *essential* or *conditionally-essential*, and quantifying the
statistical significance of these.  Although Transit is primarily
designed for Himar1 transposon libraries, some of the methods can also
be applied to other transposons, like Tn5.

**NEW in Transit 4.0:**
The GUI and command-line versions of Transit have been 
overhauled to now take combined_wig files (and corresponding metadata files),
instead of individual .wig files.
Thus, after running TPP on all the individual samples,
the multiple .wig file are combined together in a single combined_wig file
with multiple columns.  This makes it much more convenient to 
manage and analyze large datasets.


Typical Workflow
----------------

Make a flow-chart for this?...

* .fastq -> TPP -> .wig

* multiple .wig files for individual  samples -> export command -> single combined_wig file (user also writes metadata file describing samples)

* combined_wig+metadata file -> GUI or command-line Transit to run analyses -> output file (spreadsheet)

Most of the analyses in Transit produce an output file in tab-separated format that can be 
opened as a spreadsheet.  For example, 'resampling' generates an output file
that lists the data (mean insertion counts in 2 conditions being compared) and statistics (log-fold-change, p-value, adjusted p-value)
for each gene in the genome.  The user can then view this output file (in the GUI, or open it in a program like Excel)
to examine genes predicted to be conditionally essential.

(We could also show that the genome sequence (.fna or .fasta file) is needed for the TPP step,
and the genome annotation (.prot_table file) is needed for the analyses.)


TPP (Transit Pre-Processor)
---------------------------

TPP is a tool for analyzing raw sequencing data (.fastq files)
and tabulating insertion counts at TA sites.
TPP generates .wig files as output, which 
simply contain a pair of columns: 1) coordinates of TA sites,
and 2) counts of insertions observed at each site.

If Himar1 Tn libraries are prepared according to standard protocols,
such as (LONG), then the sequencing reads will have a prefix
correspond to the transposon terminus, followed by genomic DNA (along with template barcodes in read 2).
(See XXX for recommendations on sequencing parameters (length, depth), also replicates.)
TPP removes the prefix and maps the genomic portion of the read into
the genome (using BWA, [ref]).  It also reduces counts to unqiue
templates at each site, which reduces noise.

For diagnostics, each TPP run also generates a .tnseq_stats output file
which contains various metrics, such as the density (percent of TA sites
with insertions), total reads, mapped reads, and the NZmean (mean counts at non-zero sites).
It also reports statistics on how many reads were rejected because they contained
primer, adapter, or vectors sequences, which can be useful for diagnosing problems with the library.


For Tn5 (such as magellan6, see PROTOCOL), the difference in TPP are: 
1) there is no prefix, since the transposon terminus does not appear in the reads, 
and 2) the .wig files contain coordinates and counts for all sites in the 
genome, since Tn5 is not restricted to insertions at TA sites (unlike Himar1).


Combined_wig Files
------------------

Once the .wig files are created, they are combined into a combined_wig
file.  This is a new step in Transit 4.0.  Multiple .wig files are
combined into a single combined_wig file using the 'transit export
combined_wig' command on the command line.  Users then prepare a
metadata file describing each of the samples (as a spreadsheet,
e.g. in Excel, which is then saved in tab-separted format).  Most
commonly, there will be several replicates will be associated with
each condition.  The analysis methods in Transit have been redesigned
to take combined_wig and metadata files as inputs, which simplifies
things when working with large datasets.

(show command, or link to it?)

By default, the insertion counts in each dataset (.wig file) are **normalized**
by TTR (Total Trimmed Read-count) when they are combined in a combined_wig file.
This facilitates comparing insertion counts at individaul TA sites across samples.
If alternative normalization (or none) is desired, this can be specified
using a flag in the 'export combined_wig' command.

*It is often useful to examine at the pattern of insertions in conditionally-essential genes
in combined_wig files to confirm what the statistical analyses label as "significant" genes,
e.g. to ensure that the result is not biased by outliers (high counts) at individual sites,
but rather that apparent differences in counts between conditions are reflected as a consistent trend
over multiple TA sites in a gene.  This is a recommended practice.*

Once a combined_wig file is prepared, it can be used to
:ref:`assess data quality <transit_quality_control>`. There are two methods
available.  First, there is a :ref:`tnseq_stats <tnseq_stats>` command, which will
calculate various metrics for each sample, include saturation
(density, percent of TA sites with non-zero insertions), mean read
count (NZmean), as well as skewness and other statistics of the
read-count distribution (for diagnostic purposes).  This can be run at
the command-line and in the GUI. Also, plots of read-count
distributions can be generated for selected samples in the GUI (again,
helpful for identifying highly skewed samples).  A discussion about
skewed samples, the problems they cause, and what to do about them can be
found :ref:`here <transit_quality_control>`.

(also mention can do corrplots and scatterplots?)

Metadata Files
--------------

The **metadata** file describes the sample ids, filenames,
and conditions they represent (e.g. different media, growth
conditions, knockout strains, animal passaging, etc., or whatever
treatments and controls your TnSeq experiment involves).  

(metadata files are also described on the page for ZINB)

The format of the samples metadata file is a *tab-separated file* (which 
can be created/editted/saved in Excel) with 3 columns: Id, Condition, and Filename (it
must have these headers).  You can include other info about samples as additional columns, but
do not include additional rows.  Individual rows can be commented out
by prefixing them with a '#'.  Here is an example of a samples
metadata file: The filenames should match what is shown in the header
of the combined_wig (including pathnames, if present).

Note: the Condition column should have a unique label for each distinct condition (the same label shared only among replicates).
If there are attributes that distinguish the conditions (such as strain, treatment, etc), they could be included as additional columns (e.g. covariates).

Note: the filenames should match what is shown in the header of the combined_wig file;
samples are cross-referenced by filename between these two files.

::

  Id      Condition    Filename
  glyc1   glycerol     /Users/example_data/glycerol_rep1.wig
  glyc2   glycerol     /Users/example_data/glycerol_rep2.wig
  chol1   cholesterol  /Users/example_data/cholesterol_rep1.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep2.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep3.wig




Genome Annotations (prot_tables and gff files)
------------------

(this used to be documented on the Running_Transit page)

The annotation of a genome contains information about genes, such as
coordinates, strand, locus tag, gene name, and functional description.
Transit uses a custom format for annotations called "prot_table"s,
e.g. H37Rv.prot_table.  Prot_tables are **tab-separated text files**
containing the gene information in 9 specific columns:

**Prot_table file format:**

1. gene function description
2. start coordinate
3. end coordinate
4. strand
5. length of protein product (in amino acids)
6. don't care
7. don't care
8. gene name (like "dnaA")
9. ORF id (like Rv0001)

(should we put details like file formats and export/convert commands on another page???)

*It is critical that the annotation file (.prot_table) used for
analyses in Transit corresponds to exactly the same genome sequence
(.fasta or .fna) that was used to generate the .wig files with TPP,
because it is used to determine which TA sites are contained in which
genes (by coordinates).* For example, H37Rv.fna is paired with
H37Rv.prot_table, both derived from GenBank sequence NC_000962.3.

In many cases, users might often obtain annotations for their genome
in .gff (or .gff3) file format, such as downloaded from NCBI.  .gff
files contains essentially the same information about genes.  However,
there is a bit more flexibility in the .gff file format (especially in
the tags used in the right-most column), and the information about
genes is not always encoded in a uniform way, making it difficult to
use arbitrary .gffs for analyses in Transit.  Therefore, there is a
simple procedure in Transit to convert a .gff file to .prot_table
format ('**convert gff2prot_table**' via GUI or command-line).  This
step only has to be done once, and then the .prot_table can be used
for all subsequent analyses in Transit.


Command Line
------------

flags

The analysis methods in Transit are also described in this `PDF manual
<https://orca1.tamu.edu/essentiality/transit/transit-manual.pdf>`_ , focusing on 
command-line operations.




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
