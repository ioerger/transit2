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

from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import csv, misc
<<<<<<< HEAD
from pytransit.specific_tools.transit_tools import wx, basename, HAS_R, FloatVector, DataFrame, StrVector
=======
from pytransit.specific_tools.transit_tools import wx, basename
>>>>>>> a1a0f4ffbe990bbffbc1b4ac779dfcb8a82a5a95
from pytransit.globals import gui, cli, root_folder, debugging_enabled
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
            from pytransit.globals import gui, cli, root_folder, debugging_enabled
            from pytransit.specific_tools.tnseq_tools import Wig
            
            
            # 
            # get selection
            # 
            wig_objects = gui.selected_samples  or  gui.samples
            
            #
            # get read_counts and positions
            # 
            read_counts_per_wig, position_per_line = Wig.selected_as_gathered_data(wig_objects)
            number_of_wigs, number_of_lines = read_counts_per_wig.shape # => number_of_lines = len(position_per_line)
            window = 100
            for each_path_index in range(number_of_wigs):

                number_of_windows = int(number_of_lines / window) + 1  # python3 requires explicit rounding to int
                x_w = numpy.zeros(number_of_windows)
                y_w = numpy.zeros(number_of_windows)
                for window_index in range(number_of_windows):
                    x_w[window_index] = window * window_index
                    y_w[window_index] = sum(read_counts_per_wig[each_path_index][window * window_index : window * (window_index + 1)])
                
                y_smooth = stat_tools.loess(x_w, y_w, h=10000)
                plt.clf()
                plt.plot(x_w, y_w, "g+")
                plt.plot(x_w, y_smooth, "b-")
                plt.xlabel("Genomic Position (TA sites)")
                plt.ylabel("Reads per 100 insertion sites")
                
                plt.title("LOESS Fit - %s" % wig_objects[each_path_index].id)
                plt.show()