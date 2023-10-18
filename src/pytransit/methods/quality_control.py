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
    name = "Quality Control"
    
    # 
    # Quality Control
    # 
    @gui.add_wig_area_dropdown_option(name=name)
    def click_show_quality_control(event):
        from pytransit.components import qc_display
        with gui_tools.nice_error_log:
            wig_ids = [ each_sample.id for each_sample in gui.selected_samples ] 
            number_of_files = len(wig_ids)

            if number_of_files <= 0:
                raise Exception(f'''No Datasets selected, unable to run''')
            else:
                logging.log(f"Displaying results: {wig_ids}")
                try:
                    qc_window = qc_display.QualityControlFrame(gui.frame, wig_ids)
                    qc_window.Show()
                except Exception as error:
                    raise Exception(f"Error occured displaying file: {error}")