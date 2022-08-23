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

from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename, subscribe
from pytransit.universal_data import SessionData, universal
from pytransit.gui_tools import bind_to, rgba, color
from pytransit.basics.lazy_dict import LazyDict
from pytransit.components.generic.window_manager import WindowManager
from pytransit.components.generic.box import Row, Column
from pytransit.components.generic.text import Text
from pytransit.components.generic.frame import InnerFrame
from pytransit.components.annotation_area import create_annotation_area
from pytransit.components.samples_area import create_sample_area
from pytransit.components.results_area import create_results_area
from pytransit.components.parameter_panel import create_panel_area
from pytransit.components.menu import create_menu
from pytransit.analysis   import methods
from pytransit.export     import methods as export_methods
from pytransit.convert    import methods as convert_methods
from pytransit.norm_tools import methods as norm_methods

import pytransit
import pytransit.analysis
import pytransit.export
import pytransit.convert
import pytransit.components.parameter_panel as parameter_panel
import pytransit.trash as trash
import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools
import pytransit.file_display as file_display
import pytransit.qc_display as qc_display
import pytransit.images as images

class TnSeekFrame(wx.Frame):
    instructions_text = """
        1. Choose the annotation file ("prot table") that corresponds to the datasets to be analyzed.
        2. Add the desired Control and Experimental datasets.
        3. (Optional) If you wish to visualize their read counts, select the desired datasets and click on the "View" button.
        4. Select the desired analysis method from the dropdown menu on the top-right of the window, and follow its instructions.
    """.replace("\n            ","\n")
    
    # constructor
    def __init__(self, parent, DEBUG=False):
        # data accessable to all analysis methods
        universal.session_data = SessionData()
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
        
        # Timer
        self.timer = wx.Timer(self)
        def clear_status(event):
            self.statusBar.SetStatusText("")
            self.timer.Stop()
        self.Bind(wx.EVT_TIMER, clear_status, self.timer)

        
        self.SetIcon(images.transit_icon.GetIcon())

        self.workdir = os.getcwd()
        self.annotation = ""
        self.transposons = ["himar1", "tn5"]
        self.verbose = True
        
        self.statusBar = self.CreateStatusBar(1, wx.STB_SIZEGRIP, wx.ID_ANY)
        self.statusBar.SetStatusText("Welcome to TRANSIT")
        
        pub.subscribe(self.saveHistogram, "histogram")
        create_menu(self)
        
    def saveHistogram(self, msg):
        data, orf, path, delta = msg

        n, bins, patches = plt.hist(
            data, density=1, facecolor="c", alpha=0.75, bins=100
        )
        plt.xlabel("Delta Sum")
        plt.ylabel("Probability")
        plt.title("%s - Histogram of Delta Sum" % orf)
        plt.axvline(delta, color="r", linestyle="dashed", linewidth=3)
        plt.grid(True)
        genePath = os.path.join(path, orf + ".png")
        plt.savefig(genePath)
        plt.clf()

    def onHimar1Checked(self, event):
        if self.methodCheckBoxHimar1.GetValue():
            self.transposons.append("himar1")
        else:
            self.transposons.remove("himar1")
        self.filterMethodsByTransposon()

    def onTn5Checked(self, event):
        if self.methodCheckBoxTn5.GetValue():
            self.transposons.append("tn5")
        else:
            self.transposons.remove("tn5")
        self.filterMethodsByTransposon()

    def filterMethodsByTransposon(self):
        newmethods = {}
        fullmethods = pytransit.analysis.methods
        goodTn = False
        for method in fullmethods:
            goodTn = False
            for tn in self.transposons:
                if tn in fullmethods[method].transposons:
                    goodTn = True
            if goodTn:
                newmethods[method] = fullmethods[method]

        methodChoiceChoices = ["[Choose Method]"]
        for name in newmethods:
            methodChoiceChoices.append(methods[name].full_name)
        self.methodChoice.SetItems(methodChoiceChoices)
        self.methodChoice.SetSelection(0)

    def SaveFile(
        self,
        DIR=None,
        FILE="",
        WC=u'Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
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
                transit_tools.log(
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
                transit_tools.log("You chose the following file: %s" % path)
        dlg.Destroy()
        return path

    def when_loess_prev_clicked(self, event):
        from pytransit.components.samples_area import sample_table
        datasets_selected = [ each_row["path"] for each_row in sample_table.selected_rows ]
        
        if not datasets_selected:
            transit_tools.show_error_dialog("Need to select at least one control or experimental dataset.")
            return

        data, position = tnseq_tools.CombinedWig.gather_wig_data(datasets_selected)
        (K, N) = data.shape
        window = 100
        for j in range(K):

            size = (
                int(len(position) / window) + 1
            )  # python3 requires explicit rounding to int
            x_w = numpy.zeros(size)
            y_w = numpy.zeros(size)
            for i in range(size):
                x_w[i] = window * i
                y_w[i] = sum(data[j][window * i : window * (i + 1)])

            y_smooth = stat_tools.loess(x_w, y_w, h=10000)
            plt.plot(x_w, y_w, "g+")
            plt.plot(x_w, y_smooth, "b-")
            plt.xlabel("Genomic Position (TA sites)")
            plt.ylabel("Reads per 100 insertion sites")

            plt.title("LOESS Fit - %s" % transit_tools.basename(datasets_selected[j]))
            plt.show()

    def chooseNormalization(self):

        norm_methods_choices = sorted(norm_methods.keys())
        dlg = wx.SingleChoiceDialog(
            self,
            "Choose how to normalize read-counts accross datasets.",
            "Normalization Choice",
            norm_methods_choices,
            wx.CHOICEDLG_STYLE,
        )

        if dlg.ShowModal() == wx.ID_OK:
            transit_tools.log(
                "Selected the '%s' normalization method" % dlg.GetStringSelection()
            )

        dlg.Destroy()
        return dlg.GetStringSelection()
