.. transit documentation master file, created by
   sphinx-quickstart on Wed May  4 13:37:43 2016.

Welcome to TRANSIT2's documentation!
===================================

.. image:: https://img.shields.io/github/tag/ioerger/transit2.svg
    :target: https://github.com/ioerger/transit2
    :alt: GitHub last tag

Transit2 is python-based software for analyzing TnSeq data
(sequencing data from transposon mutant libraries)
to determine essentiality of bacterial genes under different conditions.

This page contains the documentation for TRANSIT2. Below are a few
quick links to some of the most important sections of the
documentation, followed by a brief overview of TRANSIT2's features.


.. NOTE::

  TRANSIT2 is a re-implementation of the `original version of TRANSIT (Transit1)
  <https://transit.readthedocs.io/en/latest/>`_,
  which is still being maintained for backwards-compatibility.
  TRANSIT2 has most of the same analytical methods, but it has an enhanced GUI.
  It relies more centrally on *combined_wig* files for storing insertion counts
  for multiple samples.
  Also, some of the command-line arguments and flags have changed.


Quick Links
~~~~~~~~~~~
.. _quick-link:


* :ref:`install-link`
* :ref:`manual-link`
* :ref:`tutorial-link`
* :ref:`tpp-link`
* :ref:`code-link`

Features
~~~~~~~~
TRANSIT2 offers a variety of features including:
 
*   More than **8 analysis methods**, including methods for determining **conditional essentiality** as well as **genetic interactions**.

*   Ability to analyze **himar1 or tn5 transposons** datasets.

*   **TrackView** to help visualize read-counts accross the genome.

*   Can **export datasets** into a variety of formats, including **IGV**.

*   Includes a **variety of normalization methods**.

*   **Quality Control** diagnostics, to idenfity poor quality datasets.

*   Ability to install as a **python package**, to import and use in your own personal scripts.   



.. 
 Mailing List
 ~~~~~~~~~~~~
 
 You can join our mailing list to get announcements of new versions, discuss any bugs, or request features! Just head over to the following site and enter your email address:
 
 
  + `https://groups.google.com/forum/#!forum/tnseq-transit/join <https://groups.google.com/forum/#!forum/tnseq-transit/join>`_
 
..


.. _manual-link:

.. toctree::
   :maxdepth: 3
   :caption: TRANSIT2 MANUAL

   transit_overview
   transit_install
   transit_running
   input_files
   quality_control

.. toctree::
   :maxdepth: 3
   :caption: [** NEW **] CRISRPi ANALYSIS METHODS 

   CGI

.. toctree::
   :maxdepth: 3
   :caption: PRE-PROCESSING

   tnseq_stats_method
   scatterplot_method
   corrplot_method
   normalize_method
   gene_means_method
   method_track_view
   method_conversions

.. toctree::
   :maxdepth: 3
   :caption: ANALYSES of single conditions

   gumbel_method
   hmm_method
   ttnfitness_method

.. toctree::
   :maxdepth: 3
   :caption: ANALYSES for pairwise comparisons

   resampling_method
   utest_method

.. toctree::
   :maxdepth: 3
   :caption: ANALYSES for multiple conditions

   gi_method
   anova_method
   zinb_method

.. toctree::
   :maxdepth: 3
   :caption: POST-PROCESSING

   pathway_enrichment_method
   heatmap_method
   method_volcano_plot


.. .. toctree::
..    :maxdepth: 3
..    :caption: Tn5 METHODS

..    method_tn5gaps


.. .. _tutorial-link:

.. .. toctree::
..    :maxdepth: 3
..    :caption: TRANSIT Tutorials

..    transit_essentiality_tutorial
..    transit_genome_tutorial
..    transit_comparative_tutorial
..    transit_interactions_tutorial
..    transit_normalization_tutorial
..    transit_export_tutorial
..    transit_console_cheatsheet

.. _tpp-link:

.. toctree::
   :maxdepth: 3
   :caption: TPP Manual

   tpp.rst

.. _code-link:

.. toctree::
   :maxdepth: 3
   :caption: Code Documentation

   transit

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Support
~~~~~~~

For any questions or comments, please contact Dr. Thomas Ioerger, ioerger@cs.tamu.edu.


