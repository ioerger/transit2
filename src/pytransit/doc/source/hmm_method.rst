.. _HMM:

HMM
===

The HMM method can be used to determine the essentiality of the entire genome, as opposed to gene-level analysis of the other methods. It is capable of identifying regions that have unusually high or unusually low read counts (i.e. growth advantage or growth defect regions), in addition to the more common categories of essential and non-essential.

.. NOTE::
   Intended only for **Himar1** datasets.

|

How does it work?
-----------------

| For a formal description of how this method works, see our paper [DeJesus2013HMM]_:
|
|  DeJesus, M.A., Ioerger, T.R. `A Hidden Markov Model for identifying essential and growth-defect regions in bacterial genomes from transposon insertion sequencing data. <http://www.ncbi.nlm.nih.gov/pubmed/24103077>`_ *BMC Bioinformatics.* 2013. 14:303


The HMM method automatically estimates the necessary internal statistical
parameters from the datasets (adjusts for global saturation 
and read counts). 


Usage
-------

::

  > python3 transit.py hmm <combined_wig_file> <metadata_file> <annotation_file> <condition_to_analyze> <output_file> [Optional Arguments]
        Optional Arguments:
            --r <string>     :=  How to handle replicates. Sum, Mean. Default: --r Mean
            --n <string>     :=  Normalization method. Default: --n TTR
            -l               :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Off.
            --iN <float>     :=  Ignore TAs occurring at given percentage (as integer) of the N terminus. Default: --iN 0
            --iC <float>     :=  Ignore TAs occurring at given percentage (as integer) of the C terminus. Default: --iC 0


Parameters
----------
You can change how the method handles
replicate datasets:

-  **Replicates:** Determines how the HMM deals with replicate datasets
   by either averaging the read-counts or summing read counts across
   datasets. For regular datasets (i.e. mean-read count > 100) the
   recommended setting is to average read-counts together. For sparse
   datasets, it summing read-counts may produce more accurate results.
|

GUI Mode
--------
The HMM analysis method can be selected from the "Method" tab in the Menu Bar. 

.. image:: _images/hmm_selection_gui.png
   :width: 1000
   :align: center

|
The parameters to input through the parameter panel for the method is equivalent to the command line usage (see parameter descriptions above for full detail): 

.. image:: _images/hmm_parameter_panel.png
   :width: 1000
   :align: center

The method is run using the combined wig, metadata, and annotation uploaded into TRANSIT.


Output and Diagnostics
----------------------

| The HMM method outputs two files. One with statistics for individual TA sites, and one
  with summaries for genes.

  The first file provides the most
  likely assignment of states for all the TA sites in the genome. Sites
  can belong to one of the following states: "E" (Essential), "GD"
  (Growth-Defect), "NE" (Non-Essential), or "GA" (Growth-Advantage). In
  addition, the output includes the probability of the particular site
  belonging to the given state. The columns of this file are defined as
  follows:

Sites Output File:
~~~~~~~~~~~~~~~~~~

+----------------+-----------------------------------------------------------------------------------------------------+
| Column Header  | Column Definition                                                                                   |
+================+=====================================================================================================+
| Location       | Coordinate of TA site                                                                               |
+----------------+-----------------------------------------------------------------------------------------------------+
| Read Count     | Observed Read Counts                                                                                |
+----------------+-----------------------------------------------------------------------------------------------------+
| Probability ES | Probability for ES state                                                                            |
+----------------+-----------------------------------------------------------------------------------------------------+
| Probability GD | Probability for GD state                                                                            |
+----------------+-----------------------------------------------------------------------------------------------------+
| Probability NE | Probability for NE state                                                                            |
+----------------+-----------------------------------------------------------------------------------------------------+
| Probability GA | Probability for GA state                                                                            |
+----------------+-----------------------------------------------------------------------------------------------------+
| State          | State Classification (ES = Essential, GD = Growth Defect, NE = Non-Essential, GA = Growth-Defect)   |
+----------------+-----------------------------------------------------------------------------------------------------+
| Gene           | Gene(s) that share(s) the TA site.                                                                  |
+----------------+-----------------------------------------------------------------------------------------------------+

|  The second file provides a gene-level classification for all the
  genes in the genome. Genes are classified as "E" (Essential), "GD"
  (Growth-Defect), "NE" (Non-Essential), or "GA" (Growth-Advantage)
  depending on the number of sites within the gene that belong to those
  states.

Genes Output File:
~~~~~~~~~~~~~~~~~~

+-------------------+-----------------------------------------------------------------------------------------------------+
| Column Header     | Column Definition                                                                                   |
+===================+=====================================================================================================+
| Orf               | Gene ID                                                                                             |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Gene Name         | Gene Name                                                                                           |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Description       | Gene Description                                                                                    |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Total Sites       | Number of TA sites                                                                                  |
+-------------------+-----------------------------------------------------------------------------------------------------+
| ES Count          | Number of sites labeled ES (Essential)                                                              |
+-------------------+-----------------------------------------------------------------------------------------------------+
| GD Count          | Number of sites labeled GD (Growth-Defect)                                                          |
+-------------------+-----------------------------------------------------------------------------------------------------+
| NE Count          | Number of sites labeled NE (Non-Essential)                                                          |
+-------------------+-----------------------------------------------------------------------------------------------------+
| GA Count          | Number of sites labeled GA (Growth-Advantage)                                                       |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Mean Insertions   | Mean insertion rate within the gene                                                                 |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Mean Reads        | Mean read count within the gene                                                                     |
+-------------------+-----------------------------------------------------------------------------------------------------+
| State Call        | State Classification (ES = Essential, GD = Growth Defect, NE = Non-Essential, GA = Growth-Defect)   |
+-------------------+-----------------------------------------------------------------------------------------------------+

|
|  Note: Libraries that are too sparse (e.g. < 30%) or which contain
  very low read-counts may be problematic for the HMM method, causing it
  to label too many Growth-Defect genes.

|

Run-time
--------

| The HMM method takes less than 10 minutes to complete. The parameters
  of the method should not affect the running-time.

|

.. rst-class:: transit_sectionend
----
