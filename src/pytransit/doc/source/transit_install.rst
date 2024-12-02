

.. _install-link:

Installation
============


Installing Transit2 using Pip
-----------------------------

You can use pip to install the TRANSIT package.

::

    > sudo pip3 install transit2

This will automatically download and install TRANSIT as a package (from PyPi), 
and all remaining required python packages. Once TRANSIT is installed as a package, 
it can be executed as a command ('transit').


.. NOTE::
   If you will be using the pre-processor, TPP, you will also need to install :ref:`install BWA <bwa-unix>`.

.. NOTE::
   The Transit package will try to install wxPython. If you encounter problems with it, see :ref:`install wxPython <install-wxpython>`

.. NOTE::
  **Note about Python3.12**: The update to Python3.12 made changes in the setuptools package that 
  disrupted the installation of Transit as well as other python packages (including some Transit
  dependencies, like wxPython).
  
  We are working on fixing these problems caused by 3.12.  In the meantime,
  if you have difficulty installing Transit with Python3.12, we recommend trying it with Python3.11, with
  which Transit should work fine.

|

Installing Transit2 using Git
-----------------------------

Alternatively, TRANSIT2 can be downloaded from the public GitHub server,
`http://github.com/ioerger/transit2 <http://github.com/ioerger/transit2>`_. 
It is released under a GPL License.  An archive 
with the lastest version of the source code can be downloaded at the following link:
`Source code.zip <https://github.com/ioerger/transit2/archive/master.zip>`_


You can clone the git respository as  follows:   

::

  > git clone https://github.com/ioerger/transit2/


TRANSIT2 is python-based. You must have **python3** installed (installed by
default on most systems). In addition, TRANSIT relies on some python 
packages/libraries/modules that you might need to install (see `Requirements`_).

You should be able to run it like this:

::

  > python3 <TRANSIT_PATH>/src/transit.py

If you encounter problems, please `contact us <https://people.engr.tamu.edu/ioerger/index.html>`_ or head to the :ref:`Troubleshooting section <transit-troubleshoot>`.


|

Requirements
------------

The following libraries/modules are required to run TRANSIT:

+ `Python 3+ <http://www.python.org>`_
+ `Numpy <http://www.numpy.org/>`_ (tested on 1.16.0)
+ `Statsmodels <https://pypi.org/project/statsmodels/>`_ (tested on 0.9.0)
+ `Scipy <http://www.scipy.org/>`_ (tested on 1.2)
+ `matplotlib <http://matplotlib.org/users/installing.html>`_ (tested on 3.0)
+ `Pillow 6.0 <https://github.com/python-pillow/Pillow>`_
+ `wxpython 4+ <http://www.wxpython.org/>`_
+ `PyPubSub 4+ <https://pypi.org/project/PyPubSub/>`_ (tested on 4.0.3)
+ `Mne <https://pypi.org/project/mne/>`_ (tested on 1.5.1)
+ `super-hash <https://pypi.org/project/super-hash/>`_ (tested on 1.4.0)
+ `ez_yaml <https://pypi.org/project/ez-yaml/>`_ (tested on 2.2.0)

All of these dependencies can be installed using the following command.

::

   > pip3 install numpy scipy pillow pypubsub matplotlib statsmodels mne wxPython super_hash ez_yaml

Pip and Python are usually preinstalled in most modern operating systems.

|

.. _install-zinb:

Additional Requirements: R (statistical analysis package)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

R is called by Transit for certain commands, such as :ref:`ZINB <zinb>`, corrplot, and heatmap.
As of now, installing R is optional, and requires these additional steps...

Additional Installation Requirements for R:

 - install `R <https://www.r-project.org/>`_ (tested on v3.5.2)
 - R packages: **MASS, pscl** (run "install.packages(MASS)" etc. in R console)
 - Python packages (for python3): rpy2 (v>=3.0) (run "pip3 install rpy2" on command line) 
 - Windows users will have to add a system variable for rpy2 to correctly work
  - To set this up, go to **Control Panel** > **System** > **Advanced system settings** > **Environment Variables**, then add a new variable named ``R_HOME`` with the path to your R installation (e.g. ``C:\Program Files\R\R-4.0.5``).
 
 .. - Python packages (for python2.7): rpy2 (v<2.9.0) (run "pip install 'rpy2<2.9.0' " on command line)



.. Use as a Python Package
.. -----------------------------------------------------


.. TRANSIT can be (optionally) installed as a python package. This can simplify the installation process as it will automatically install most of the requirements. In addition, it will allow users to use some of transit functions in their own scripts if they desire. Below is a brief example of importing transit functions into python. In this example, pair of .wig files are parsed into their read-counts (data) and genomic positions (position), and then normalization factors are calculated. See the documentation of the package for further examples:

.. 

..         >>> import pytransit.norm_tools as norm_tools
..         >>> import pytransit.tnseq_tools as tnseq_tools
..         >>> (data, position) = tnseq_tools.get_data(["transit/data/cholesterol_glycerol.transit/glycerol_rep1.wig", "transit/data/cholesterol_glycerol.transit/glycerol_rep2.wig"])
..         >>> print(data)
..         array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
..                [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
..         >>> factors = norm_tools.TTR_factors(data)
..         >>> print(factors)
..         array([[ 1.        ],
..                [ 0.62862886]])



Optional: Install BWA to use with TPP pre-processor
---------------------------------------------------

If you will be using the pre-processor, TPP, you will also need to install `BWA <http://bio-bwa.sourceforge.net/>`_.




.. _bwa-unix:

Linux & OSX Instructions
~~~~~~~~~~~~~~~~~~~~~~~~

Download the source files:


 + `http://sourceforge.net/projects/bio-bwa/files/ <http://sourceforge.net/projects/bio-bwa/files/>`_


Extract the files:

::


    > tar -xvjf bwa-0.7.12.tar.bz2


Go to the directory with the extracted source-code, and run make to create the executable files:

::


    > cd bwa-0.7.12
    > make


.. _bwa-win:

Windows Instructions
~~~~~~~~~~~~~~~~~~~~

For Windows, we provide a windows executable (.exe) for Windows 64 bit:

  + `bwa-0.7.12_windows.zip <http://saclab.tamu.edu/essentiality/transit/bwa-0.7.12_windows.zip>`_



The 32-bit version of Windows is not recommended as it is limited in the amount of system memory that can be used.


|

.. _transit-upgrade:

Upgrading
---------

The process of upgrading transit will depend on how you installed transit initially.


Method 1: Upgrading package installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


If you installed TRANSIT as a package, then to upgrade, simply use pip to install tnseq-transit again, but this time include the '--upgrade' flag. For example:


::

    > sudo pip install transit2 --upgrade

This will automatically download and install the latest version of TRANSIT, as well as upgrade any of its requirements if necessary for compatability.


Method 2: Upgrading source installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you installed TRANSIT by downloading the raw source, then you can upgrade TRANSIT simply by replacing the old source code with the latest version. You can obtain a .zip archive with the latest version of the source through the following link:

https://github.com/ioerger/transit2/archive/master.zip

Simply extract the code, and replace your existing files or delete the directory with the old source doe and use the newest version.

Or you may perform a 'git pull'

|

.. NOTE::
   If an an older version of wxPython is already installed (< 4.0), you may have to remove it and install version 4.0+.

|

.. _install-wxpython:

Installing wxPython
-------------------

wxPython 4+ can be manually installed using pip

::

   > pip3 install wxPython

If the above command fails and you already have wxPython < 4.0 installed, you may have to manually remove it.
See https://stackoverflow.com/questions/50688630/cannot-uninstall-wxpython-3-0-2-0-macos for details.

.. NOTE::

  Installing *wxPython* can be a bit finicky.  It might require installing the
  development version of GTK first.  There are at least two versions currently, 
  *gtk2* and *gtk3*.
  Transit should work with both, although there can be small differences in the 
  visual look of the GUI.  To get *wxPython* to install, you might try doing this:

    > sudo apt-get install libgtk-2-dev

    or

    > sudo apt-get install libgtk-3-dev

  depending on which version of *libgtk* you have installed.

.. NOTE::

  If you are still having problems, another solution might be to
  install the appropriate version of wxPython directly from the
  developer's website.  For example, on a recent Linux machine, I had
  to manually do this step *before* doing 'pip install transit2' (in a
  fresh virtual environment):

    > pip install -U -f https://extras.wxpython.org/wxPython4/extras/linux/gtk3/ubuntu-20.04 wxPython

  But you might have to adjust the URL to the version specific for your environment.

.. _transit-troubleshoot:

Troubleshooting
---------------

1. No window appears when running in GUI mode.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


This problem is likely due to running OSX and previously unsuported versions of matplotlib.
Please upgrade matplotlib to the latest version using:

::

    > pip3 install 'matplotlib' --upgrade

|

2. pip3: SystemError: Cannot compile 'Python.h'.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs when you do not have the development libraries for python. You can fix this by installing the python-dev packages:


::

    > sudo apt-get install python-dev


|

3. pip: "The following required packages can not be built: freetype,png," etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs when you do not have some dependencies that are necessary to build some of the python modules TRANSIT requires (usually matplotlib). Installing the following linux dependencies should fix this:

::

    > sudo apt-get install libpng-dev libjpeg8-dev libfreetype6-dev


|

4. pip3: "No lapack/blas resources found"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs when you do not have some dependencies that are necessary to build some of the python modules TRANSIT requires (usually numpy/scipy). Installing the following linux dependencies should fix this:


::

    > sudo apt-get install libblas-dev liblapack-dev libatlas-base-dev gfortran


|

5. "resources.ContextualVersionConflict (six 1.5.2)..."
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs some of the python modules are out of date. You can use pip to upgrade them as follows:


::

    > sudo pip3 install six --upgrade
