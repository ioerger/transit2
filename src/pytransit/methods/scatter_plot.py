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
from pytransit.specific_tools.transit_tools import wx, basename, HAS_R, FloatVector, DataFrame, StrVector, EOL
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components import samples_area, file_display, results_area, parameter_panel
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    name = "Scatter Plot"
    
    # TODO: confirm menu option not needed
    # @gui.add_menu("Preprocessing", name)
    # def on_menu_click(event):
    #     pass
    
    # 
    # Scatter Plot
    # 
    @samples_area.create_sample_area_button(name=name, size=(120, -1))
    def create_show_scatter_plot_button(event):
        with gui_tools.nice_error_log:
            import numpy
            import matplotlib
            import matplotlib.pyplot as plt
            from pytransit.specific_tools import stat_tools
            selected_samples = gui.selected_samples
            if len(selected_samples) == 2:
                logging.log( f"Showing scatter plot for: {[ each_sample.id for each_sample in selected_samples ]}")
                from pytransit.specific_tools.transit_tools import gather_sample_data_for
                data, position = gather_sample_data_for(selected_samples=True)
                x = data[0, :]
                y = data[1, :]

                plt.plot(x, y, "bo")
                plt.title("Scatter plot - Reads at TA sites")
                plt.xlabel(selected_samples[0].id)
                plt.ylabel(selected_samples[1].id)
                plt.show()
            else:
                # NOTE: was a popup
                logging.error(f"Select 2 samples (not {len(selected_samples)})")