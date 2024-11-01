# TRANSIT

[![Version](https://img.shields.io/github/tag/ioerger/transit2.svg)](https://github.com/ioerger/transit2)   [![Build Status](https://travis-ci.org/ioerger/transit2.svg?branch=master)](https://travis-ci.org/ioerger/transit2)   [![Documentation Status](https://readthedocs.org/projects/transit2/badge/?version=latest)](http://transit2.readthedocs.io/en/latest/?badge=latest)   [![Downloads](https://pepy.tech/badge/transit2)](https://pepy.tech/project/transit2)


**NOTE: TRANSIT2 requires python3.6+.**

Welcome! This is the distribution for the TRANSIT2 and TPP tools developed by the [Ioerger Lab](http://orca2.tamu.edu/tom/iLab.html) at Texas A&M University.

TRANSIT2 is a tool for processing and statistical analysis of Tn-Seq data.
It provides an easy to use graphical interface and access to three different analysis methods that allow the user to determine essentiality in a single condition as well as between conditions.

TRANSIT2 Home page: http://saclab.tamu.edu/essentiality/transit/index_transit2.html

TRANSIT2 Documentation: https://transit2.readthedocs.io/en/stable/transit_overview.html

[Changelog](https://github.com/ioerger/transit2/blob/master/CHANGELOG.md)


## Features
TRANSIT2 offers a variety of features including:

-   More than **10 analysis methods**, including methods for determining **conditional essentiality** as well as **genetic interactions**.

-   Ability to analyze datasets from libraries constructed using  **himar1 or tn5 transposons**.

-   **TrackView** to help visualize read-counts across the genome.

-   Can **export datasets** into a variety of formats, including **IGV**.

-   Includes a **variety of normalization methods**.

-   **Quality Control** diagnostics, to idenfity poor quality datasets.

-   Ability to install as a **python package**, to import and use in your own personal scripts.



## Support

For any questions or comments, please contact Dr. Thomas Ioerger, ioerger@cs.tamu.edu.




## Instructions

For full instructions on how to install and run TRANSIT2 (and the optional pre-processor, TPP), 
please see the documentation included in this distribution ("src/pytransit/doc" folder) or visit the following web page:


https://transit2.readthedocs.io/en/stable/


## Datasets

The TRANSIT distribution comes with some example .wig files in the data/ directory, as well as an example annotation file (.prot\_table format) in the data/genomes/ directory. Additional genomes may be found on the following website:

http://saclab.tamu.edu/essentiality/transit/genomes/


## Copyright Information

Source code for TRANSIT2 and TPP are available open source under the terms of the GNU General Public License (Version 3.0) as published by the Free Software Foundation. For more information on this license, please see the included LICENSE.md file or visit their website at:

http://www.gnu.org/licenses/gpl.html

<!-- NOTE: sudo apt install -y libeigen3-dev -->
<!-- NOTE: sudo apt install -y 'python-wxgtk3*' -->
<!-- pip install --upgrade setuptools -->
<!-- pip3 install attrdict -->