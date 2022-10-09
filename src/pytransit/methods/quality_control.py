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
import pytransit.generic_tools.csv as csv
import pytransit.generic_tools.misc as misc
from pytransit.specific_tools.transit_tools import wx, basename, HAS_R, FloatVector, DataFrame, StrVector
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components import samples_area, file_display, results_area, parameter_panel
from pytransit.components.spreadsheet import SpreadSheet

@misc.singleton
class Method:
    name = "Quality Control"
    
    # TODO: confirm menu option not needed
    # @gui.add_menu("Preprocessing", name)
    # def on_menu_click(event):
    #     pass
    
    # 
    # Quality Control
    # 
    @samples_area.create_sample_area_button(name=name, size=(130,-1))
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