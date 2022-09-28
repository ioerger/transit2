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

from pytransit.methods.analysis                  import methods
from pytransit.methods.convert                   import methods as convert_methods
from pytransit.methods.export                    import methods as export_methods
from pytransit.tools.gui_tools                   import bind_to, rgba, color
from pytransit.tools.norm_tools                  import methods as norm_methods
from pytransit.tools.transit_tools               import HAS_WX, wx, GenBitmapTextButton, pub, basename, subscribe
from pytransit.components.generic.box            import Row, Column
from pytransit.components.generic.frame          import InnerFrame
from pytransit.components.generic.text           import Text
from pytransit.components.generic.window_manager import WindowManager
from pytransit.components.menu                   import create_menu
from pytransit.components.parameter_panel        import create_panel_area
from pytransit.components.results_area           import create_results_area
from pytransit.components.samples_area           import create_sample_area
from pytransit.components.annotation_area        import create_annotation_area
from pytransit.basics.lazy_dict                  import LazyDict
from pytransit.universal_data                    import SessionData, universal

from pytransit.tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, stat_tools
import pytransit
import pytransit.components.parameter_panel as parameter_panel
import pytransit.components.trash as trash
import pytransit.components.file_display as file_display
import pytransit.components.qc_display as qc_display
import pytransit.components.images as images

class TnSeqFrame(wx.Frame):
    instructions_text = """
        1. Choose the annotation file ("prot table") that corresponds to the datasets to be analyzed.
        2. Add the desired Control and Experimental datasets.
        3. (Optional) If you wish to visualize their read counts, select the desired datasets and click on the "View" button.
        4. Select the desired analysis method from the dropdown menu on the top-right of the window, and follow its instructions.
    """.replace("\n            ","\n")
    
    # constructor
    def __init__(self, parent, DEBUG=False):
        # data accessable to all analysis methods
        universal.frame = self
        # connect to GUI tools (otherwise they will not function)
        gui_tools.bit_map = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_OTHER, (16, 16))
        
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
                    
                    # children
                    if True:
                        data_column.add(
                            create_annotation_area(self),
                            proportion=0,
                        )
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

        self.workdir = os.getcwd()
        self.transposons = ["himar1", "tn5"]
        self.verbose = True
        
        self.status_bar = self.CreateStatusBar(1, wx.STB_SIZEGRIP, wx.ID_ANY)
        self.status_bar.SetStatusText("Welcome to TRANSIT")
        
        pub.subscribe(self.save_histogram, "histogram")
        create_menu(self)
        
    def save_histogram(self, msg):
        data, orf, path, delta = msg
        n, bins, patches = plt.hist(data, density=1, facecolor="c", alpha=0.75, bins=100)
        plt.xlabel("Delta Sum")
        plt.ylabel("Probability")
        plt.title("%s - Histogram of Delta Sum" % orf)
        plt.axvline(delta, color="r", linestyle="dashed", linewidth=3)
        plt.grid(True)
        genePath = os.path.join(path, orf + ".png")
        plt.savefig(genePath)
        plt.clf()

    def SaveFile(
        self,
        DIR=None,
        FILE="",
        WC='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
    ):
        """
        Create and show the Save FileDialog
        """
        path = ""

        if not DIR:
            DIR = os.getcwd()

        
        dlg = wx.FileDialog(
            self,
            message="Save file as ...",
            defaultDir=DIR,
            defaultFile=FILE,
            wildcard=WC,
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
        )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if self.verbose:
                logging.log(
                    "You chose the following output filename: %s" % path
                )
        dlg.Destroy()
        return path

    def OpenFile(self, DIR=".", FILE="", WC=""):
        """
        Create and show the Open FileDialog
        """
        path = ""
        
        dlg = wx.FileDialog(
            self,
            message="Save file as ...",
            defaultDir=DIR,
            defaultFile=FILE,
            wildcard=WC,
            style=wx.FD_OPEN,
        )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if self.verbose:
                logging.log("You chose the following file: %s" % path)
        dlg.Destroy()
        return path

    def choose_normalization(self):
        norm_methods_choices = sorted(norm_methods.keys())
        dlg = wx.SingleChoiceDialog(
            self,
            "Choose how to normalize read-counts accross datasets.",
            "Normalization Choice",
            norm_methods_choices,
            wx.CHOICEDLG_STYLE,
        )

        if dlg.ShowModal() == wx.ID_OK:
            logging.log(
                "Selected the '%s' normalization method" % dlg.GetStringSelection()
            )

        dlg.Destroy()
        return dlg.GetStringSelection()
