

.. _tnseq_stats:

.. rst-class:: transit_clionly
tnseq_stats
===========

You can generate the same table to statistics as on the Quality Control panel in the GUI
from the command-line using the 'tnseq_stats' command.  Here is an example:

::

  > python3 src/transit.py tnseq_stats --help

  usage: python3 src/transit.py tnseq_stats <file.wig>+ [-o <output_file>]
         python3 src/transit.py tnseq_stats -c <combined_wig> [-o <output_file>]

  > python3 src/transit.py tnseq_stats -c src/pytransit/data/cholesterol_glycerol.transit/comwig.tsv

  dataset density mean_ct NZmean  NZmedian        max_ct  total_cts       skewness        kurtosis
  src/pytransit/data/cholesterol_glycerol.transit/cholesterol_rep1.wig   0.44    139.6   317.6   147     125355.5        10414005.0      54.8    4237.7
  src/pytransit/data/cholesterol_glycerol.transit/cholesterol_rep2.wig   0.44    171.4   390.5   148     704662.8        12786637.9      105.8   14216.2
  src/pytransit/data/cholesterol_glycerol.transit/cholesterol_rep3.wig   0.36    173.8   484.2   171     292294.8        12968502.500000002      42.2    2328.0
  src/pytransit/data/cholesterol_glycerol.transit/glycerol_rep1.wig      0.42    123.3   294.5   160     8813.3  9195672.4       4.0     33.0
  src/pytransit/data/cholesterol_glycerol.transit/glycerol_rep2.wig      0.52    123.8   240.1   127     8542.5  9235984.2       4.0     33.5


.. rst-class:: transit_sectionend
----
