# -*- coding: utf-8 -*-
# Copyright 2015.
#   Michael A. DeJesus, Chaitra Ambadipudi, and  Thomas R. Ioerger.
#
#
#    This file is part of TRANSIT.
#
#    TRANSIT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License.
#
#
#    TRANSIT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with TRANSIT.  If not, see <http://www.gnu.org/licenses/>.

from collections import defaultdict
from functools import partial
import datetime
import math
import multiprocessing as mp
import os
import subprocess
import sys
import threading
import time
import traceback

import numpy
import matplotlib
import matplotlib.pyplot as plt

from pytransit.specific_tools.gui_tools                   import bind_to, rgba, color
from pytransit.specific_tools.norm_tools                  import methods as norm_methods
from pytransit.specific_tools.transit_tools               import HAS_WX, wx, GenBitmapTextButton, basename
from pytransit.components.generic.box            import Row, Column
from pytransit.components.generic.frame          import InnerFrame
from pytransit.components.generic.text           import Text
from pytransit.components.generic.window_manager import WindowManager
from pytransit.components.menu                   import create_menu
from pytransit.components.parameter_panel        import create_panel_area
from pytransit.components.results_area           import create_results_area
from pytransit.components.samples_area           import create_sample_area
from pytransit.generic_tools.lazy_dict                  import LazyDict
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

from pytransit.specific_tools import gui_tools, transit_tools, tnseq_tools, norm_tools, stat_tools
import pytransit
import pytransit.components.parameter_panel as parameter_panel
import pytransit.components.trash as trash
import pytransit.components.file_display as file_display
import pytransit.components.images as images

class TnSeqFrame(wx.Frame):
    # constructor
    def __init__(self, parent):
        # data accessable to all analysis methods
        gui.frame = self
        # connect to GUI tools (otherwise they will not function)
        gui_tools.bit_map = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_OTHER, (16, 16))
        
        gui.debug_wx_if_needed()
        
        with InnerFrame(parent, title="TRANSIT") as frame:
            
            with Row() as main_wrapper:
                
                with Column() as spacer:
                    spacer.add(Text(" "))
                    main_wrapper.add(
                        spacer,
                        proportion=0,
                    )
                    
                # 
                # data column
                # 
                with Column() as data_column:
                    data_column.wx_object.Add(10,5) # padding
                    # children
                    if True:
                        data_column.add(
                            create_sample_area(self),
                            proportion=1,
                            expand=True,
                            border=5,
                        )
                        data_column.add(
                            create_results_area(self),
                            expand=True,
                            border=5,
                        )
                    
                    main_wrapper.add(
                        data_column,
                        expand=True,
                        proportion=5,
                        border=5,
                    )
                
                # 
                # panel column
                # 
                with Column() as panel_column:
                    
                    if True:
                        panel_column.add(
                            create_panel_area(self),
                            expand=True,
                            border=5,
                        )
                    
                    main_wrapper.add(
                        panel_column,
                        expand=True,
                        proportion=2,
                    )
            
                frame.add(
                    main_wrapper,
                    expand=True,
                )
            
            self.inner_frame = frame

        self.Centre(wx.BOTH)
        self.SetIcon(images.transit_icon.GetIcon())
        
        self.status_bar = self.CreateStatusBar(1, wx.STB_SIZEGRIP, wx.ID_ANY)
        self.status_bar.SetStatusText("Welcome to TRANSIT")
        
        create_menu(self)

# 
# PNG render options
# 
@transit_tools.ResultsFile
class PngFile:
    @staticmethod
    def can_load(path):
        return path.endswith(".png")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        from pytransit.specific_tools import gui_tools
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type="Image",
            path=self.path,
            # could potentially add file date here
            
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Image": lambda *args: gui_tools.show_image(self.path),
            })
        )