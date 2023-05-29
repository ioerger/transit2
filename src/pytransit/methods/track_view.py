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
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    name = "Track View"
    
    # 
    # Track View
    # 
    @staticmethod
    @gui.add_wig_area_dropdown_option(name=name)
    def click_show_track_view(event):
        with gui_tools.nice_error_log:
            import pytransit.components.trash as trash
            annotation_path = gui.annotation_path
            wig_ids = [ each_sample.id for each_sample in gui.selected_samples ]

            if wig_ids and annotation_path:
                if debugging_enabled:
                    logging.log("Visualizing counts for: %s" % ", ".join(wig_ids))
                view_window = trash.TrashFrame(gui.frame, wig_ids, annotation_path, gene="")
                view_window.Show()
            elif not wig_ids:
                # NOTE: was a popup
                logging.error("Error: No samples selected.")
            else:
                # NOTE: was a popup
                logging.error("Error: No annotation file selected.")