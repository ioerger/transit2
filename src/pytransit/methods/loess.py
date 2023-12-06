import sys
import os
import time
import ntpath
import math
import random
import datetime
import collections
import heapq

import numpy

from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import csv, misc
from pytransit.specific_tools.transit_tools import wx, basename
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

@misc.singleton
class Method:
    name = "LOESS Plot"
    
    # 
    # LOESS
    # 
    @gui.add_wig_area_dropdown_option(name=name)
    def click_show_loess(event):
        with gui_tools.nice_error_log:
            import numpy
            import matplotlib
            import matplotlib.pyplot as plt
            from pytransit.specific_tools import stat_tools
            from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
            from pytransit.specific_tools.tnseq_tools import Wig
            from statsmodels.nonparametric.smoothers_lowess import lowess
            
            # 
            # get selection
            # 
            wig_objects = gui.selected_samples  or  gui.samples
            
            #
            # get read_counts and positions
            # 
            read_counts_per_wig, position_per_line = Wig.selected_as_gathered_data(wig_objects)
            TAsites = position_per_line 
            counts = numpy.mean(read_counts_per_wig,axis=0)
            win = 1000 # 1kb
            rounded_ta_sites = (TAsites/win).astype(int)*win
            buckets = {}
            for b, cnt in zip(rounded_ta_sites, counts):
                if b not in buckets:
                    buckets[b] = []
                buckets[b].append(cnt)
            sorted_keys = tuple(sorted(buckets.keys()))
            Y = [numpy.mean(buckets[key]) for key in sorted_keys]
            X = tuple(sorted_keys)
            Y_normalized = stat_tools.loess_correction(X=X,Y=Y)

            fit = lowess(exog=X,endog=Y_normalized) # using statsmodels; note: previously, stat_tools.loess(x_w, y_w, h=10000) was being called
            fit_no_correction = lowess(exog=X,endog=Y)

            fig, ((ax1, ax3)) = plt.subplots(1, 2)
            # fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2)
            fig.set_figheight(5)
            fig.set_figwidth(10)

            ax1.plot(X,Y, "g+")
            ax1.plot(fit_no_correction[:,0],fit_no_correction[:,1], "b-")
            ax1.set_xlabel("Genomic Position (bp)")
            ax1.set_ylabel("(log-scale) mean insertion count over %sbp windows" % win)
            ax1.set_yscale("log")
            ax1.set_title("Before Correction")

            # ax2.plot(X,Y, "g+")
            # ax2.plot(fit_no_correction[:,0],fit_no_correction[:,1], "b-")
            # ax2.set_xlabel("Genomic Position (bp)")
            # ax2.set_ylabel("mean insertion count over %sbp windows" % win)
            # ax2.set_title("Fit Before Correction")
            
            ax3.plot(X,Y, "y+")
            ax3.plot(fit[:,0],fit[:,1], "b-")
            ax3.set_xlabel("Genomic Position (bp)")
            # ax3.set_ylabel("mean insertion count over %sbp windows" % win)
            ax3.set_yscale("log")
            ax3.set_title("After LOESS Correction")

            # ax4.plot(X,Y, "y+")
            # ax4.plot(fit[:,0],fit[:,1], "b-")
            # ax4.set_xlabel("Genomic Position (bp)")
            # ax4.set_ylabel("mean insertion count over %sbp windows" % win)
            # ax4.set_title("LOESS Fit")
            
            plt.show()

