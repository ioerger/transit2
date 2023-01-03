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
from pytransit.components import samples_area, file_display, results_area, parameter_panel
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    name = "LOESS"
    
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
            buckets = {}
            for i in range(len(TAsites)):
              coord,cnt = TAsites[i],counts[i]
              b = int(coord/float(win))
              if b not in buckets: buckets[b] = []
              buckets[b].append(cnt)
            X = sorted(buckets.keys())
            Y = [numpy.mean(buckets[x]) for x in X]
            Xbp = [x*win for x in X]
            fit = lowess(exog=Xbp,endog=Y) # using statsmodels; note: previously, stat_tools.loess(x_w, y_w, h=10000) was being called

            fig, (ax1, ax2) = plt.subplots(1, 2)
            fig.set_figheight(5)
            fig.set_figwidth(10)

            ax1.plot(Xbp,Y, "g+")
            ax1.plot(fit[:,0],fit[:,1], "b-")
            ax1.set_xlabel("Genomic Position (bp)")
            ax1.set_ylabel("mean insertion count over %sbp windows" % win)
            ax1.set_yscale("log")
            ax1.set_title("LOESS Fit (log-scale)")

            ax2.plot(Xbp,Y, "g+")
            ax2.plot(fit[:,0],fit[:,1], "b-")
            ax2.set_xlabel("Genomic Position (bp)")
            ax2.set_ylabel("mean insertion count over %sbp windows" % win)
            ax2.set_title("LOESS Fit")
            m,s = numpy.max(fit[:,1]),numpy.std(Y)
            ax2.set_ylim(0,m+2*s)

            plt.show()

