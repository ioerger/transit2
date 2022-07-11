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
import ez_yaml

from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename, subscribe
from pytransit.core_data import SessionData, universal
from pytransit.gui_tools import bind_to, rgba, color
from pytransit.basics.lazy_dict import LazyDict
from pytransit.components.comwig_picker import create_comwig_picker
from pytransit.components.generic.table import Table
import pytransit
import pytransit.analysis
import pytransit.export
import pytransit.convert
import pytransit.trash as trash
import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools
import pytransit.file_display as file_display
import pytransit.qc_display as qc_display
import pytransit.images as images

methods           = pytransit.analysis.methods
export_methods    = pytransit.export.methods
convert_methods   = pytransit.convert.methods
norm_methods      = norm_tools.methods

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
        
        # connect to GUI tools (otherwise they will not function)
        gui_tools.window  = self
        gui_tools.bit_map = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_OTHER, (16, 16))
        
        
        wx.Frame.__init__(
            self,
            parent,
            id=wx.ID_ANY,
            title=u"TRANSIT",
            pos=wx.DefaultPosition,
            size=wx.Size(1350, 975),
            style=wx.DEFAULT_FRAME_STYLE | wx.TAB_TRAVERSAL,
        )

        # 
        # main window
        # 
        if True:
            self.mainWindow = wx.ScrolledWindow(
                self,
                wx.ID_ANY,
                wx.DefaultPosition,
                wx.Size(-1, -1),
                wx.HSCROLL | wx.VSCROLL,
            )
            self.mainWindow.SetScrollRate(5, 5)
            self.mainWindow.SetMinSize(wx.Size(700, -1))
            
            # 
            # windowWrapper
            # 
            if True:
                windowWrapper = wx.BoxSizer(wx.HORIZONTAL)
                # 
                # windowSizer
                # 
                # if True:
                #     windowSizer = wx.BoxSizer(wx.VERTICAL)
                    
                #     # 
                #     # orgSizer
                #     # 
                #     if True:
                #         orgSizer = wx.StaticBoxSizer(
                #             wx.StaticBox(
                #                 self.mainWindow,
                #                 wx.ID_ANY,
                #                 u"Organism"
                #             ),
                #             wx.VERTICAL,
                #         )
                        
                #         # 
                #         # annotation sizer
                #         # 
                #         if True:
                #             annot_sizer = wx.BoxSizer(wx.HORIZONTAL)
                            
                #             # 
                #             # text
                #             # 
                #             if True:
                #                 label_annot = wx.StaticText(
                #                     self.mainWindow,
                #                     wx.ID_ANY,
                #                     u"Annotation File:",
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     0,
                #                 )
                #                 annot_sizer.Add(
                #                     label_annot,
                #                     0,
                #                     wx.ALIGN_CENTER_VERTICAL,
                #                     0
                #                 )
                            
                #             # 
                #             # picker
                #             # 
                #             if True:
                #                 self.annotationFilePicker = wx.FilePickerCtrl(
                #                     self.mainWindow,
                #                     id=wx.ID_ANY,
                #                     size=(400, 30),
                #                     wildcard=u"prot_table or GFF3 files (*.gff3;*.gff;*.prot_table;*.txt)|*.gff3;*.gff;*.prot_table;*.txt",
                #                     message="Select Annotation file (.prot_table or .gff3)",
                #                     style=wx.FLP_DEFAULT_STYLE | wx.FLP_USE_TEXTCTRL | wx.FD_MULTIPLE,
                #                 )
                #                 self.annotationFilePicker.SetInitialDirectory(os.getcwd())
                #                 annot_sizer.Add(
                #                     self.annotationFilePicker,
                #                     proportion=1,
                #                     flag=wx.EXPAND | wx.ALL,
                #                     border=5,
                #                 )
                #                 self.annotationFilePicker
                                
                #                 @bind_to(self.annotationFilePicker, wx.EVT_FILEPICKER_CHANGED)
                #                 def annotationFileFunc(self, event):
                #                     self.annotation = event.GetPath()

                #             orgSizer.Add(annot_sizer, 1, wx.EXPAND, 5)

                #         windowSizer.Add(orgSizer, 0, wx.EXPAND, 5)

                #     # 
                #     # Samples
                #     # 
                #     if True:
                #         ctrlSizer = wx.StaticBoxSizer(
                #             wx.StaticBox(
                #                 self.mainWindow,
                #                 wx.ID_ANY,
                #                 u"Samples"
                #             ),
                #             wx.VERTICAL,
                #         )
                        
                #         # 
                #         # box
                #         # 
                #         if True:
                #             ctrlBoxSizer2 = wx.BoxSizer(wx.HORIZONTAL)
                            
                #             # 
                #             # ctrlViewButton
                #             # 
                #             if True:
                #                 self.ctrlViewButton = wx.Button(
                #                     self.mainWindow,
                #                     wx.ID_ANY,
                #                     u"Track View",
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     0,
                #                 )
                #                 self.ctrlViewButton.Hide()
                #                 ctrlBoxSizer2.Add(self.ctrlViewButton, 0, wx.ALL, 5)
                                
                #                 @bind_to(self.ctrlViewButton, wx.EVT_BUTTON)
                #                 def _(event):
                #                     self.allViewFunc(event)
                                

                #             # 
                #             # combined_wig_file_picker
                #             # 
                #             if True:
                #                 self.combined_wig_file_picker = create_comwig_picker(window=self)
                #                 ctrlBoxSizer2.Add(self.combined_wig_file_picker, 1, wx.ALIGN_CENTER_VERTICAL, 5)
                                
                #             ctrlSizer.Add(ctrlBoxSizer2, 0, wx.EXPAND, 5)
                        
                #         # 
                #         # wig_table
                #         # 
                #         with Table() as (wx_object, component):
                #             self.wig_table = wx_object
                #             ctrlSizer.Add(wx_object, 1, wx.ALL | wx.EXPAND, 5)

                #         windowSizer.Add(
                #             ctrlSizer,
                #             1,
                #             wx.EXPAND,
                #             5
                #         )
                    
                #     # 
                #     # conditionsSizer
                #     # 
                #     if True:
                #         conditionsSizer = wx.StaticBoxSizer(
                #             wx.StaticBox(self.mainWindow, wx.ID_ANY, u"Conditions"),
                #             wx.VERTICAL,
                #         )
                        
                #         # 
                #         # box
                #         # 
                #         if True:
                #             boxSizer = wx.BoxSizer(wx.HORIZONTAL)
                            
                #             # 
                #             # experimentTrackViewButton
                #             # 
                #             if True:
                #                 self.experimentTrackViewButton = wx.Button(
                #                     self.mainWindow,
                #                     wx.ID_ANY,
                #                     u"Track View",
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     0,
                #                 )
                #                 self.experimentTrackViewButton.Hide()
                #                 boxSizer.Add(self.experimentTrackViewButton, 0, wx.ALL, 5)
                                
                #                 @bind_to(self.experimentTrackViewButton, wx.EVT_BUTTON)
                #                 def _(event):
                #                     return self.allViewFunc(event)

                #             conditionsSizer.Add(boxSizer, 0, wx.EXPAND, 5)
                        
                #         # 
                #         # Conditions
                #         # 
                #         with Table() as (wx_object, component):
                #             wx_object.SetMaxSize(wx.Size(-1, 200))
                #             conditionsSizer.Add(wx_object, 1, wx.ALL | wx.EXPAND, 5)
                #             self.conditions_table = obj
                        
                #         windowSizer.Add(conditionsSizer, 1, wx.EXPAND, 5)

                #     # 
                #     # Results
                #     # 
                #     if True:
                #         results_sizer = wx.StaticBoxSizer(
                #             wx.StaticBox(self.mainWindow, wx.ID_ANY, u"Results Files"), wx.VERTICAL
                #         )
                        
                #         # 
                #         # Box
                #         # 
                #         if True:
                #             boxSizer4 = wx.BoxSizer(wx.HORIZONTAL)
                            
                #             # 
                #             # displayButton
                #             # 
                #             if True:
                #                 self.displayButton = wx.Button(
                #                     self.mainWindow,
                #                     wx.ID_ANY,
                #                     u"Display Table",
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     0,
                #                 )
                #                 boxSizer4.Add(self.displayButton, 0, wx.ALL, 5)
                                
                #                 @bind_to(self.displayButton, wx.EVT_BUTTON)
                #                 def displayFileFunc(event):
                #                     next = self.results_table.GetNextSelected(-1)
                #                     if next > -1:
                #                         dataset = self.results_table.GetItem(next, 3).GetText()
                #                         if self.verbose:
                #                             transit_tools.transit_message(
                #                                 "Displaying results: %s"
                #                                 % self.results_table.GetItem(next, 0).GetText()
                #                             )

                #                         try:
                #                             fileWindow = file_display.TransitGridFrame(self, dataset)
                #                             fileWindow.Show()
                #                         except Exception as e:
                #                             transit_tools.transit_message(
                #                                 "Error occurred displaying file: %s" % str(e)
                #                             )
                #                             traceback.print_exc()

                #                     else:
                #                         if self.verbose:
                #                             transit_tools.transit_message("No results selected to display!")

                                
                #             # 
                #             # fileActionButton
                #             # 
                #             if True:
                #                 self.fileActionButton = wx.Button(
                #                     self.mainWindow,
                #                     wx.ID_ANY,
                #                     u"Display Graph",
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     0,
                #                 )
                #                 self.fileActionButton.Hide()

                #                 boxSizer4.Add(self.fileActionButton, 0, wx.ALL, 5)
                                
                #                 @bind_to(self.fileActionButton, wx.EVT_BUTTON)
                #                 def _(event):
                #                     self.fileActionFunc(event)
                                
                                
                                
                #             # 
                #             # addFileButton
                #             # 
                #             if True:
                #                 self.addFileButton = GenBitmapTextButton(
                #                     self.mainWindow,
                #                     1,
                #                     bit_map,
                #                     "Add Results File",
                #                     size=wx.Size(150, 30)
                #                 )
                #                 boxSizer4.Add(self.addFileButton, 0, wx.ALL, 5)
                                
                #                 @subscribe("file")
                #                 def add_file(data):
                #                     fullpath = data["path"]
                #                     type     = data["type"]
                #                     date     = data["date"]
                #                     name = transit_tools.basename(fullpath)
                #                     # FIXME
                #                     self.results_table.InsertItem(self.index_file, name)
                #                     self.results_table.SetItem(self.index_file, 1, f"{type}")
                #                     self.results_table.SetItem(self.index_file, 2, f"{date}")
                #                     self.results_table.SetItem(self.index_file, 3, f"{fullpath}")
                #                     self.index_file += 1
                                
                #                 @bind_to(self.addFileButton, wx.EVT_BUTTON)
                #                 def _(event):
                #                     try:
                #                         dlg = wx.FileDialog(
                #                             self,
                #                             message="Choose a file",
                #                             defaultDir=self.workdir,
                #                             defaultFile="",
                #                             wildcard=u"Results Files (*.dat)|*.dat;|\nResults Files (*.txt)|*.txt;|\nAll files (*.*)|*.*",
                #                             style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                #                         )
                #                         if dlg.ShowModal() == wx.ID_OK:
                #                             paths = dlg.GetPaths()
                #                             print("You chose the following Results file(s):")
                #                             for fullpath in paths:
                #                                 print("\t%s" % fullpath)
                #                                 name = transit_tools.basename(fullpath)
                #                                 line = open(fullpath).readline()
                #                                 if line.startswith("#Gumbel"):
                #                                     type = "Gumbel"
                #                                 elif line.startswith("#Binomial"):
                #                                     type = "Binomial"
                #                                 elif line.startswith("#HMM - Sites"):
                #                                     type = "HMM - Sites"
                #                                 elif line.startswith("#HMM - Genes"):
                #                                     type = "HMM - Genes"
                #                                 elif line.startswith("#Resampling"):
                #                                     type = "Resampling"
                #                                 elif line.startswith("#DE-HMM - Sites"):
                #                                     type = "DE-HMM - Sites"
                #                                 elif line.startswith("#DE-HMM - Segments"):
                #                                     type = "DE-HMM - Segments"
                #                                 elif line.startswith("#GI"):
                #                                     type = "GI"
                #                                 else:
                #                                     type = "Unknown"
                #                                 data = {
                #                                     "path": fullpath,
                #                                     "type": type,
                #                                     "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p"),
                #                                 }
                #                                 wx.CallAfter(pub.sendMessage, "file", data=data)
                #                         dlg.Destroy()
                #                     except Exception as e:
                #                         transit_tools.transit_message("Error: %s" % e)
                #                         print("PATH", fullpath)
                #                         exc_type, exc_obj, exc_tb = sys.exc_info()
                #                         fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                #                         print(exc_type, fname, exc_tb.tb_lineno)
                                
                #             # 
                #             # fileActionChoice
                #             # 
                #             if True:
                #                 fileActionChoiceChoices = [  u"[Choose Action]"  ]
                #                 self.fileActionChoice = wx.Choice(
                #                     self.mainWindow,
                #                     wx.ID_ANY,
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     fileActionChoiceChoices,
                #                     0,
                #                 )
                #                 self.fileActionChoice.SetSelection(0)
                #                 boxSizer4.Add(self.fileActionChoice, 0, wx.ALL, 5)
                                
                #                 @bind_to(self.fileActionChoice, wx.EVT_CHOICE)
                #                 def fileActionFunc(event):
                #                     # 0 - nothing
                #                     # 1 - Volcano
                #                     # 2 - Hist gene counts ratio
                #                     plot_choice = self.fileActionChoice.GetCurrentSelection()
                #                     plot_name = self.fileActionChoice.GetString(plot_choice)
                #                     if plot_name == "[Choose Action]":
                #                         return
                #                     next = self.results_table.GetNextSelected(-1)
                #                     if next > -1:
                #                         dataset_path = self.results_table.GetItem(next, 3).GetText()
                #                         dataset_name = self.results_table.GetItem(next, 0).GetText()
                #                         dataset_type = self.results_table.GetItem(next, 1).GetText()

                #                         if self.verbose:
                #                             transit_tools.transit_message(
                #                                 "Performing the '%s' action on dataset '%s'"
                #                                 % (plot_name, dataset_name)
                #                             )

                #                         if plot_name == "Create a Volcano Plot":
                #                             self.graphVolcanoPlot(dataset_name, dataset_type, dataset_path)
                #                         elif plot_name == "Plot Histogram of logFC of Gene Counts":
                #                             self.graphGeneCounts(dataset_name, dataset_type, dataset_path)
                #                         elif plot_name == "Plot Ranked Probability of Essentiality":
                #                             self.graphRankedZbar(dataset_name, dataset_type, dataset_path)
                #                         else:
                #                             return

                #                         self.fileActionChoice.SetSelection(0)
                #                     else:
                #                         transit_tools.ShowError(MSG="Please select a results file to plot!")

                                
                #             results_sizer.Add(boxSizer4, 0, 0, 5)
                        
                #         # 
                #         # results_table
                #         # 
                #         with Table(
                #             initial_columns=[ "Name", "Type", "Date", "Full Path"],
                #             max_size=(-1, 200)
                #         ) as (wx_object, component):
                            
                #             self.results_table = obj
                #             results_sizer.Add(wx_object, 1, wx.ALL | wx.EXPAND, 5)
                            
                #     windowSizer.Add(results_sizer, 1, wx.EXPAND, 5)
                #     self.mainWindow.SetSizer(windowSizer)
                #     self.mainWindow.Layout()
                #     windowSizer.Fit(self.mainWindow)
                #     windowWrapper.Add(self.mainWindow, 1, wx.ALL | wx.EXPAND, 5)
                
                # # 
                # # m_panel5
                # # 
                # if True:
                #     self.m_panel5 = wx.Panel(
                #         self,
                #         wx.ID_ANY,
                #         wx.DefaultPosition,
                #         wx.DefaultSize,
                #         wx.DOUBLE_BORDER | wx.TAB_TRAVERSAL,
                #     )
                #     self.m_panel5.SetMaxSize(wx.Size(2, -1))
                #     windowWrapper.Add(self.m_panel5, 0, wx.ALL | wx.EXPAND, 5)
                
                # # 
                # # options window
                # # 
                # if True:
                #     self.optionsWindow = wx.ScrolledWindow(
                #         self,
                #         wx.ID_ANY,
                #         wx.DefaultPosition,
                #         wx.Size(-1, -1),
                #         wx.HSCROLL | wx.VSCROLL | wx.EXPAND,
                #     )
                #     self.optionsWindow.SetScrollRate(5, 5)
                #     self.optionsWindow.SetMinSize(wx.Size(310, 1000))
                    
                #     # 
                #     # box
                #     # 
                #     if True:
                #         self.optionsSizer = wx.BoxSizer(wx.VERTICAL)
                        
                #         # 
                #         # Logo Section
                #         # 
                #         if True:
                #             self.logoImg = wx.StaticBitmap(
                #                 self.optionsWindow,
                #                 wx.ID_ANY,
                #                 wx.NullBitmap,
                #                 wx.DefaultPosition,
                #                 wx.DefaultSize,
                #                 0,
                #             )
                #             self.optionsSizer.Add(self.logoImg, 0, wx.ALL | wx.ALIGN_CENTER, 5)
                        
                #         # 
                #         # versionLabel
                #         # 
                #         if True:
                #             self.versionLabel = wx.StaticText(
                #                 self.optionsWindow,
                #                 wx.ID_ANY,
                #                 u"",
                #                 wx.DefaultPosition,
                #                 (100, 25),
                #                 wx.ALIGN_CENTRE,
                #             )
                #             self.versionLabel.Wrap(-1)
                #             self.versionLabel.SetFont(wx.Font(10, 74, 90, 92, False, "Sans"))

                #             self.optionsSizer.Add(
                #                 self.versionLabel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                #             )
                        
                #         # 
                #         # methodInfoSizer
                #         # 
                #         if True:
                #             self.methodInfoText = wx.StaticBox(self.optionsWindow, wx.ID_ANY, u"Instructions")
                #             self.methodInfoText.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
                #             self.methodInfoSizer = wx.StaticBoxSizer(self.methodInfoText, wx.VERTICAL)
                            
                #             # 
                #             # methodShortText
                #             # 
                #             if True:
                #                 self.methodShortText = wx.StaticText(
                #                     self.optionsWindow, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
                #                 )
                #                 self.methodShortText.Wrap(250)
                #                 self.methodShortText.Hide()
                #                 self.methodInfoSizer.Add(
                #                     self.methodShortText, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                #                 )
                            
                #             # 
                #             # methodLongText
                #             # 
                #             if True:
                #                 self.methodLongText = wx.StaticText(
                #                     self.optionsWindow, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
                #                 )
                #                 self.methodLongText.Wrap(250)
                #                 self.methodLongText.Hide()
                #                 self.methodInfoSizer.Add(
                #                     self.methodLongText, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                #                 )
                            
                #             # 
                #             # methodDescText
                #             # 
                #             if True:

                #                 self.methodDescText = wx.StaticText(
                #                     self.optionsWindow, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
                #                 )
                #                 self.methodDescText.Wrap(250)
                #                 self.methodDescText.Hide()
                #                 self.methodInfoSizer.Add(
                #                     self.methodDescText, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                #                 )
                            
                #             # 
                #             # methodTnText
                #             # 
                #             if True:
                #                 self.methodTnText = wx.StaticText(
                #                     self.optionsWindow, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
                #                 )
                #                 self.methodTnText.Wrap(250)

                #                 font = wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.BOLD)
                #                 self.methodTnText.SetFont(font)
                #                 self.methodTnText.Hide()

                #                 self.methodInfoSizer.Add(
                #                     self.methodTnText, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                #                 )
                            
                #             # 
                #             # methodInstructions
                #             # 
                #             if True:
                #                 self.methodInstructions = wx.StaticText(
                #                     self.optionsWindow,
                #                     wx.ID_ANY,
                #                     self.instructions_text,
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     0,
                #                 )
                #                 self.methodInstructions.Wrap(250)
                #                 self.methodInfoSizer.Add(
                #                     self.methodInstructions, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                #                 )
                            
                #             self.optionsSizer.Add(self.methodInfoSizer, 0, wx.ALL | wx.EXPAND, 5)
                        
                #         # 
                #         # Method Options
                #         # 
                #         if True:
                #             self.methodSizerText = wx.StaticBox(self.optionsWindow, wx.ID_ANY, u"Method Options")
                #             self.methodSizerText.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
                #             self.methodSizer = wx.StaticBoxSizer(self.methodSizerText, wx.VERTICAL)
                            
                #             # 
                #             # methodPanel1
                #             # 
                #             if True:
                #                 self.methodPanel1 = wx.Panel(
                #                     self.optionsWindow,
                #                     wx.ID_ANY,
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     wx.TAB_TRAVERSAL,
                #                 )
                #                 self.methodPanel1.SetMinSize(wx.Size(50, 1))
                #                 self.methodSizer.Add(self.methodPanel1, 0, wx.ALL, 5)
                            
                #             # 
                #             # globalLabel
                #             # 
                #             if True:
                #                 self.globalLabel = wx.StaticText(
                #                     self.optionsWindow,
                #                     wx.ID_ANY,
                #                     u"Global Options",
                #                     wx.DefaultPosition,
                #                     (130, 20),
                #                     0,
                #                 )
                #                 self.globalLabel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
                #                 self.methodSizer.Add(self.globalLabel, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)
                            
                #             # 
                #             # globalPanel
                #             # 
                #             if True:
                #                 self.globalPanel = wx.Panel(
                #                     self.optionsWindow,
                #                     wx.ID_ANY,
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     wx.TAB_TRAVERSAL,
                #                 )
                #                 self.globalPanel.SetMinSize(wx.Size(280, 150))
                #                 self.globalPanel.SetMaxSize(wx.Size(-1, -1))

                #                 globalSizerVT = wx.BoxSizer(wx.VERTICAL)
                #                 nTermSizer    = wx.BoxSizer(wx.HORIZONTAL)
                #                 cTermSizer    = wx.BoxSizer(wx.HORIZONTAL)

                #                 # N TERMINUS - GLOBAL
                #                 self.globalNTerminusLabel = wx.StaticText(
                #                     self.globalPanel,
                #                     wx.ID_ANY,
                #                     u"Ignore N-Terminus %:",
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     0,
                #                 )
                #                 self.globalNTerminusLabel.Wrap(-1)
                #                 self.globalNTerminusText = wx.TextCtrl(self.globalPanel, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0)
                #                 self.globalNTerminusIcon = pytransit.analysis.base.InfoIcon(
                #                     self.globalPanel,
                #                     wx.ID_ANY,
                #                     tooltip="Ignores a fraction of the ORF, beginning at the N-terminal end. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",
                #                 )
                #                 nTermSizer.Add(self.globalNTerminusLabel, 1, wx.ALIGN_CENTER, 5)
                #                 nTermSizer.Add(self.globalNTerminusText , 1, wx.ALIGN_CENTER, 5)
                #                 nTermSizer.Add(self.globalNTerminusIcon , 1, wx.ALIGN_CENTER, 5)

                #                 # C TERMINUS - GLOBAL
                #                 self.globalCTerminusLabel = wx.StaticText(
                #                     self.globalPanel,
                #                     wx.ID_ANY,
                #                     u"Ignore C-Terminus %:",
                #                     wx.DefaultPosition,
                #                     wx.DefaultSize,
                #                     0,
                #                 )
                #                 self.globalCTerminusLabel.Wrap(-1)
                #                 self.globalCTerminusText = wx.TextCtrl(self.globalPanel, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0)
                #                 self.globalCTerminusIcon = pytransit.analysis.base.InfoIcon(
                #                     self.globalPanel,
                #                     wx.ID_ANY,
                #                     tooltip="Ignores a fraction of the ORF, beginning at the C-terminal end. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",
                #                 )

                #                 cTermSizer.Add(self.globalCTerminusLabel, 1, wx.ALIGN_CENTER_VERTICAL, 5)
                #                 cTermSizer.Add(self.globalCTerminusText, 1, wx.ALIGN_CENTER_VERTICAL, 5)
                #                 cTermSizer.Add(self.globalCTerminusIcon, 1, wx.ALIGN_CENTER, 5)

                #                 # Control Libraries text - GLOBAL
                #                 ctrlLibSizer = wx.BoxSizer(wx.HORIZONTAL)
                #                 self.ctrlLibLabel = wx.StaticText(
                #                     self.globalPanel,
                #                     wx.ID_ANY,
                #                     u"Control Libraries:",
                #                     wx.DefaultPosition,
                #                     (170, -1),
                #                     0,
                #                 )
                #                 self.ctrlLibLabel.Wrap(-1)
                #                 self.ctrlLibText = wx.TextCtrl(self.globalPanel, wx.ID_ANY, "", wx.DefaultPosition, (-1, -1), 0)
                #                 self.ctrlLibTip = pytransit.analysis.base.InfoIcon(
                #                     self.globalPanel,
                #                     wx.ID_ANY,
                #                     tooltip="String of letters representing an \
                #                     identifier for the libraries the datasets belong to. For example, if adding three \
                #                     datasets of different libraries, change the string to 'ABC'. Set of letters used  \
                #                     must match those in Experimental datasets. Keep empty or with all letters equal, e.g. \
                #                     'AAA', to do regular resampling.",
                #                 )

                #                 self.ctrlLibText.Disable()
                #                 ctrlLibSizer.Add(self.ctrlLibLabel, 0, wx.ALIGN_CENTER_VERTICAL, 5)
                #                 ctrlLibSizer.Add(self.ctrlLibText, 0, wx.ALIGN_CENTER_VERTICAL, 5)
                #                 ctrlLibSizer.Add(self.ctrlLibTip, 0, wx.ALIGN_CENTER_VERTICAL, 5)

                #                 # Experimental Libraries text - GLOBAL
                #                 expLibSizer = wx.BoxSizer(wx.HORIZONTAL)
                #                 self.expLibLabel = wx.StaticText(
                #                     self.globalPanel,
                #                     wx.ID_ANY,
                #                     u"Experimental Libraries:",
                #                     wx.DefaultPosition,
                #                     (170, -1),
                #                     0,
                #                 )
                #                 self.expLibLabel.Wrap(-1)
                #                 self.expLibText = wx.TextCtrl(
                #                     self.globalPanel, wx.ID_ANY, "", wx.DefaultPosition, (-1, -1), 0
                #                 )
                #                 self.expLibTip = pytransit.analysis.base.InfoIcon(
                #                     self.globalPanel,
                #                     wx.ID_ANY,
                #                     tooltip="String of letters representing an identifier for the libraries the datasets \
                #                     belong to. For example, if adding three datasets of different libraries, change the \
                #                     string to 'ABC'. Set  of letters used must match those in Control datasets. Keep \
                #                     empty or with all letters equal, e.g. 'AAA', to do regular resampling.",
                #                 )

                #                 self.expLibText.Disable()
                #                 expLibSizer.Add(self.expLibLabel, 0, wx.ALIGN_CENTER_VERTICAL, 5)
                #                 expLibSizer.Add(self.expLibText , 0, wx.ALIGN_CENTER_VERTICAL, 5)
                #                 expLibSizer.Add(self.expLibTip  , 0, wx.ALIGN_CENTER_VERTICAL, 5)

                #                 globalSizerVT.Add(nTermSizer  , 1, wx.EXPAND, 5)
                #                 globalSizerVT.Add(cTermSizer  , 1, wx.EXPAND, 5)
                #                 globalSizerVT.Add(ctrlLibSizer, 1, wx.EXPAND, 5)
                #                 globalSizerVT.Add(expLibSizer , 1, wx.EXPAND, 5)

                #                 self.globalPanel.SetSizer(globalSizerVT)
                #                 self.globalPanel.Layout()
                #                 globalSizerVT.Fit(self.globalPanel)
                #                 self.methodSizer.Add(self.globalPanel, 1, wx.ALIGN_CENTER_HORIZONTAL, 5)
                            
                #         self.optionsSizer.Add(self.methodSizer, 0, wx.EXPAND, 5)

                #     self.optionsWindow.SetSizer(self.optionsSizer)
                #     self.optionsWindow.Layout()

                #     self.optionsWindow.Fit()

                # windowWrapper.Add(self.optionsWindow, 0, wx.ALL, 5)
        # --------------------#
        
        # TODO: cleanup formatting of nested items
        self.SetSizer(windowWrapper)
        self.Layout()
        # if True:
        #     self.m_menubar1 = wx.MenuBar(0)
        #     self.fileMenuItem = wx.Menu()
        #     self.exportMenuItem = wx.Menu()
        #     self.selectedExportMenuItem = wx.Menu()

        #     # Selected datasets
        #     self.exportMenuItem.AppendSubMenu(
        #         self.selectedExportMenuItem, u"Selected Datasets"
        #     )

        #     self.fileMenuItem.AppendSubMenu(self.exportMenuItem, u"Export")

        #     self.convertMenuItem = wx.Menu()
        #     self.annotationConvertPTToPTTMenu = wx.MenuItem(
        #         self.convertMenuItem,
        #         wx.ID_ANY,
        #         u"prot_table to PTT",
        #         wx.EmptyString,
        #         wx.ITEM_NORMAL,
        #     )
        #     self.convertMenuItem.Append(self.annotationConvertPTToPTTMenu)

        #     self.annotationConvertPTToGFF3Menu = wx.MenuItem(
        #         self.convertMenuItem,
        #         wx.ID_ANY,
        #         u"prot_table to GFF3",
        #         wx.EmptyString,
        #         wx.ITEM_NORMAL,
        #     )
        #     self.convertMenuItem.Append(self.annotationConvertPTToGFF3Menu)

        #     self.annotationConvertPTTToPT = wx.MenuItem(
        #         self.convertMenuItem,
        #         wx.ID_ANY,
        #         u"PTT to prot_table",
        #         wx.EmptyString,
        #         wx.ITEM_NORMAL,
        #     )

        #     self.convertMenuItem.Append(self.annotationConvertPTTToPT)

        #     # self.annotationConvertGFF3ToPT = wx.MenuItem( self.convertMenuItem, wx.ID_ANY, u"GFF3 to prot_table", wx.EmptyString, wx.ITEM_NORMAL )
        #     # self.convertMenuItem.Append( self.annotationConvertGFF3ToPT )
        #     self.fileMenuItem.AppendSubMenu(self.convertMenuItem, u"Convert")

        #     self.fileExitMenuItem = wx.MenuItem(
        #         self.fileMenuItem, wx.ID_ANY, u"&Exit", wx.EmptyString, wx.ITEM_NORMAL
        #     )
        #     self.fileMenuItem.Append(self.fileExitMenuItem)
        #     self.m_menubar1.Append(self.fileMenuItem, u"&File")

        #     self.viewMenuItem = wx.Menu()
        #     self.scatterMenuItem = wx.MenuItem(
        #         self.viewMenuItem,
        #         wx.ID_ANY,
        #         u"&Scatter Plot",
        #         wx.EmptyString,
        #         wx.ITEM_NORMAL,
        #     )

        #     self.viewMenuItem.Append(self.scatterMenuItem)

        #     self.trackMenuItem = wx.MenuItem(
        #         self.viewMenuItem, wx.ID_ANY, u"&Track View", wx.EmptyString, wx.ITEM_NORMAL
        #     )

        #     self.viewMenuItem.Append(self.trackMenuItem)
        #     self.m_menubar1.Append(self.viewMenuItem, u"&View")

        #     #
        #     self.methodsMenuItem = wx.Menu()
        #     self.himar1MenuItem = wx.Menu()
        #     self.tn5MenuItem = wx.Menu()

        #     self.methodsMenuItem.AppendSubMenu(self.himar1MenuItem, "&Himar1 Methods")
        #     self.methodsMenuItem.AppendSubMenu(self.tn5MenuItem, "&Tn5 Methods")
        #     self.m_menubar1.Append(self.methodsMenuItem, u"&Analysis")

        #     self.SetMenuBar(self.m_menubar1)

        #     self.helpMenuItem = wx.Menu()
        #     self.documentationMenuItem = wx.MenuItem(
        #         self.helpMenuItem,
        #         wx.ID_ANY,
        #         u"&Documentation",
        #         wx.EmptyString,
        #         wx.ITEM_NORMAL,
        #     )
        #     self.helpMenuItem.Append(self.documentationMenuItem)
        #     self.aboutMenuItem = wx.MenuItem(
        #         self.helpMenuItem, wx.ID_ANY, u"&About", wx.EmptyString, wx.ITEM_NORMAL
        #     )
        #     self.helpMenuItem.Append(self.aboutMenuItem)

        #     self.m_menubar1.Append(self.helpMenuItem, u"&Help")

        #     self.statusBar = self.CreateStatusBar(1, wx.STB_SIZEGRIP, wx.ID_ANY)

        self.Centre(wx.BOTH)
        
        # 
        # Connect Menu Events
        # 
        # if True:
        #     self.Bind(wx.EVT_MENU, self.annotationPT_to_PTT , id=self.annotationConvertPTToPTTMenu.GetId(),  )
        #     self.Bind(wx.EVT_MENU, self.annotationPT_to_GFF3, id=self.annotationConvertPTToGFF3Menu.GetId(), )
        #     self.Bind(wx.EVT_MENU, self.annotationPTT_to_PT , id=self.annotationConvertPTTToPT.GetId(),      )
        #     self.Bind(wx.EVT_MENU, self.Exit                , id=self.fileExitMenuItem.GetId()               )
        #     self.Bind(wx.EVT_MENU, self.scatterFunc         , id=self.scatterMenuItem.GetId()                )
        #     self.Bind(wx.EVT_MENU, self.allViewFunc         , id=self.trackMenuItem.GetId()                  )
        #     self.Bind(wx.EVT_MENU, self.aboutFunc           , id=self.aboutMenuItem.GetId()                  )
        #     self.Bind(wx.EVT_MENU, self.documentationFunc   , id=self.documentationMenuItem.GetId()          )
            
        #     self.timer = wx.Timer(self)
        #     self.Bind(wx.EVT_TIMER, self.clearStatus, self.timer)

        
        # self.SetIcon(images.transit_icon.GetIcon())

        # self.workdir = os.getcwd()
        # self.annotation = ""
        # self.transposons = ["himar1", "tn5"]

        # self.logoImg.SetBitmap(images.transit_logo2.GetImage().ConvertToBitmap())
        # self.versionLabel.SetLabel(pytransit.__version__)
        # self.methodSizerText.Hide()


        # self.index_file = 0
        

        # self.verbose = True

        # self.statusBar.SetStatusText("Welcome to TRANSIT")
        # self.progress_count = 0
        # pub.subscribe(self.setProgressRange, "progressrange")
        # pub.subscribe(self.updateProgress, "progress")
        # pub.subscribe(self.updateStatus, "status")
        # pub.subscribe(self.finishRun, "finish")
        # pub.subscribe(self.saveHistogram, "histogram")

        # # 
        # # Export Menu Items
        # # 
        # for name in export_methods:
        #     export_methods[name].gui.defineMenuItem(self, export_methods[name].label)
        #     tempMenuItem = export_methods[name].gui.menuitem
        #     self.selectedExportMenuItem.Append(tempMenuItem)

        #     self.Bind(
        #         wx.EVT_MENU,
        #         partial(self.ExportSelectFunc, export_methods[name].label),
        #         tempMenuItem,
        #     )

        # # Convert Menu Items
        # for name in convert_methods:
        #     convert_methods[name].gui.defineMenuItem(self, convert_methods[name].label)
        #     tempMenuItem = convert_methods[name].gui.menuitem
        #     self.convertMenuItem.Append(tempMenuItem)

        #     self.Bind(
        #         wx.EVT_MENU,
        #         partial(self.ConvertSelectFunc, convert_methods[name].label),
        #         tempMenuItem,
        #     )

        # # Method Panels

        # methodChoiceChoices = ["[Choose Method]"]
        # methodorder = [("gumbel", 1), ("resampling", 2), ("hmm", 3)]
        # order = defaultdict(lambda: 100)
        # for k, v in methodorder:
        #     order[k] = v

        # for name in sorted(methods.keys(), key=lambda x: order[x]):
        #     methods[name].gui.definePanel(self)
        #     # methods[name].gui.panel.BackgroundColour = (0, 200, 20)
        #     self.methodSizer.Add(
        #         methods[name].gui.panel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
        #     )
        #     methods[name].gui.Hide()

        #     if "himar1" in methods[name].transposons:
        #         tempMenuItem = wx.MenuItem(
        #             self.himar1MenuItem,
        #             wx.ID_ANY,
        #             methods[name].fullname(),
        #             wx.EmptyString,
        #             wx.ITEM_NORMAL,
        #         )
        #         self.Bind(
        #             wx.EVT_MENU,
        #             partial(self.MethodSelectFunc, methods[name].fullname()),
        #             tempMenuItem,
        #         )

        #         self.himar1MenuItem.Append(tempMenuItem)

        #     if "tn5" in methods[name].transposons:
        #         tempMenuItem = wx.MenuItem(
        #             self.tn5MenuItem,
        #             wx.ID_ANY,
        #             methods[name].fullname(),
        #             wx.EmptyString,
        #             wx.ITEM_NORMAL,
        #         )
        #         self.Bind(
        #             wx.EVT_MENU,
        #             partial(self.MethodSelectFunc, methods[name].fullname()),
        #             tempMenuItem,
        #         )
        #         self.tn5MenuItem.Append(tempMenuItem)

        # # progress
        # self.progressPanel = wx.Panel(
        #     self.optionsWindow,
        #     wx.ID_ANY,
        #     wx.DefaultPosition,
        #     wx.DefaultSize,
        #     wx.TAB_TRAVERSAL,
        # )
        # progressSizer = wx.BoxSizer(wx.VERTICAL)

        # self.progressLabel = wx.StaticText(
        #     self.progressPanel,
        #     wx.ID_ANY,
        #     u"Progress",
        #     wx.DefaultPosition,
        #     wx.DefaultSize,
        #     0,
        # )
        # self.progressLabel.Wrap(-1)
        # progressSizer.Add(self.progressLabel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        # self.progress = wx.Gauge(
        #     self.progressPanel,
        #     wx.ID_ANY,
        #     20,
        #     wx.DefaultPosition,
        #     wx.DefaultSize,
        #     wx.GA_HORIZONTAL | wx.GA_SMOOTH,
        # )
        # progressSizer.Add(self.progress, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        # # self.progressPanel.BackgroundColour = (0, 0, 250)
        # self.progressPanel.SetSizer(progressSizer)
        # self.progressPanel.SetMaxSize(wx.Size(100, 100))
        # self.progressPanel.Layout()

        # # progressSizer.Fit( self.progressPanel )
        # self.methodSizer.Add(
        #     self.progressPanel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
        # )

        # # self.methodSizer.Add( self.globalLabel, 1, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        # # self.methodSizer.Hide()
        # self.progress.SetRange(50)
        # #########

        # self.optionsWindow.Fit()

        # self.HideProgressSection()
        # self.HideGlobalOptions()

    def Exit(self, event):
        """Exit Menu Item"""
        if self.verbose:
            transit_tools.transit_message("Exiting Transit")
        self.Close()

    def updateProgress(self, msg):
        """"""
        method, count = msg
        self.progress_count = count
        try:
            self.progress.SetValue(self.progress_count)
        except:
            pass

    def setProgressRange(self, msg):
        """"""
        count = msg
        try:
            self.progress.SetRange(count)
        except:
            pass

    def updateStatus(self, msg, time=-1):
        """"""
        if type(msg) == type("A"):
            text = msg
        else:
            method, text, time = msg
        if time > 0:
            self.timer.Start(time)
        self.statusBar.SetStatusText(text)

    def clearStatus(self, event):
        self.statusBar.SetStatusText("")
        self.timer.Stop()

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

    def finishRun(self, msg):
        try:
            self.progress_count = 0
            self.progress.SetValue(self.progress_count)

        except Exception as e:
            transit_tools.transit_message("Error: %s" % e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def ResetProgress(self):
        self.progress_count = 0

    def HideAllOptions(self):
        self.HideGlobalOptions()
        self.HideProgressSection()
        for name in methods:
            methods[name].gui.Hide()

    def HideGlobalOptions(self):
        self.globalLabel.Hide()
        self.HideGlobalCTerminus()
        self.HideGlobalNTerminus()
        self.HideGlobalLibraries()

    def HideGlobalCTerminus(self):
        self.globalCTerminusLabel.Hide()
        self.globalCTerminusText.Hide()
        self.globalCTerminusIcon.Hide()

    def HideGlobalNTerminus(self):
        self.globalNTerminusLabel.Hide()
        self.globalNTerminusText.Hide()
        self.globalNTerminusIcon.Hide()

    def HideGlobalLibraries(self):
        self.HideGlobalCtrlLibraries()
        self.HideGlobalExpLibraries()

    def HideGlobalCtrlLibraries(self):
        self.ctrlLibLabel.Hide()
        self.ctrlLibText.Hide()
        self.ctrlLibTip.Hide()

    def HideGlobalExpLibraries(self):
        self.expLibLabel.Hide()
        self.expLibText.Hide()
        self.expLibTip.Hide()

    def ShowGlobalOptions(self):
        self.globalLabel.Show()
        self.ShowGlobalCTerminus()
        self.ShowGlobalNTerminus()
        self.ShowGlobalLibraries()

    def ShowGlobalCTerminus(self):
        self.globalCTerminusLabel.Show()
        self.globalCTerminusText.Show()
        self.globalCTerminusIcon.Show()

    def ShowGlobalNTerminus(self):
        self.globalNTerminusText.Show()
        self.globalNTerminusLabel.Show()
        self.globalNTerminusIcon.Show()

    def ShowGlobalLibraries(self):
        self.ShowGlobalCtrlLibraries()
        self.ShowGlobalExpLibraries()

    def ShowGlobalCtrlLibraries(self):
        self.ctrlLibLabel.Show()
        self.ctrlLibText.Show()
        self.ctrlLibTip.Show()

    def ShowGlobalExpLibraries(self):
        self.expLibLabel.Show()
        self.expLibText.Show()
        self.expLibTip.Show()

    def HideProgressSection(self):
        self.progressLabel.Hide()
        self.progress.Hide()

    def ShowProgressSection(self):
        self.progressLabel.Show()
        self.progress.Show()

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
            methodChoiceChoices.append(methods[name].fullname())
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
                transit_tools.transit_message(
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
                transit_tools.transit_message("You chose the following file: %s" % path)
        dlg.Destroy()
        return path

    
    def allSelected(self, col=5):
        selected_all = self.ctrlSelected(col) + self.expSelected(col)
        return selected_all

    def loadCtrlFile(self, fullpath):
        name = transit_tools.basename(fullpath)
        (
            density,
            meanrd,
            nzmeanrd,
            nzmedianrd,
            maxrd,
            totalrd,
            skew,
            kurtosis,
        ) = tnseq_tools.get_wig_stats(fullpath)
        self.wig_table.InsertItem(self.index_ctrl, name)
        self.wig_table.SetItem(self.index_ctrl, 1, "%1.1f" % (totalrd))
        self.wig_table.SetItem(self.index_ctrl, 2, "%2.1f" % (density * 100))
        self.wig_table.SetItem(self.index_ctrl, 3, "%1.1f" % (meanrd))
        self.wig_table.SetItem(self.index_ctrl, 4, "%d" % (maxrd))
        self.wig_table.SetItem(self.index_ctrl, 5, "%s" % (fullpath))

        self.wig_table.Select(self.index_ctrl)
        self.index_ctrl += 1
        try:
            self.ctrlLibText.SetValue(self.ctrlLibText.GetValue() + "A")
        except Exception as e:
            transit_tools.transit_message("Error Modifying Ctrl Lib String: %s" % e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def ctrlRemoveFunc(self, event):
        next = self.wig_table.GetNextSelected(-1)
        while next != -1:
            if self.verbose:
                transit_tools.transit_message(
                    "Removing control item (%d): %s"
                    % (next, self.wig_table.GetItem(next, 0).GetText())
                )

            # Update library string after removing wig file
            updated_lib_text = self.ctrlLibText.GetValue()
            updated_lib_text = updated_lib_text[:next] + updated_lib_text[(next + 1) :]
            self.ctrlLibText.SetValue(updated_lib_text)

            # Delete and Get next selected
            self.wig_table.DeleteItem(next)
            next = self.wig_table.GetNextSelected(-1)
            self.index_ctrl -= 1

    def allViewFunc(self, event, gene=""):

        annotationpath = self.annotation
        datasets = self.ctrlSelected() + self.expSelected()

        if datasets and annotationpath:
            if self.verbose:
                transit_tools.transit_message(
                    "Visualizing counts for: %s"
                    % ", ".join([transit_tools.fetch_name(d) for d in datasets])
                )
            viewWindow = trash.TrashFrame(self, datasets, annotationpath, gene=gene)
            viewWindow.Show()
        elif not datasets:
            transit_tools.ShowError("Error: No datasets selected.")
            return
        else:
            transit_tools.ShowError("Error: No annotation file selected.")
            return

    def scatterFunc(self, event):
        """ """
        # annotationpath = self.annotation
        datasets = self.ctrlSelected() + self.expSelected()
        if len(datasets) == 2:
            if self.verbose:
                transit_tools.transit_message(
                    "Showing scatter plot for: %s"
                    % ", ".join([transit_tools.fetch_name(d) for d in datasets])
                )
            (data, position) = tnseq_tools.get_data(datasets)
            X = data[0, :]
            Y = data[1, :]

            plt.plot(X, Y, "bo")
            plt.title("Scatter plot - Reads at TA sites")
            plt.xlabel(transit_tools.fetch_name(datasets[0]))
            plt.ylabel(transit_tools.fetch_name(datasets[1]))
            plt.show()
        else:
            transit_tools.ShowError(
                MSG="Please make sure only two datasets are selected (across control and experimental datasets)."
            )

    def aboutFunc(self, event):
        description = """TRANSIT is a tool for analysing TnSeq data. It provides an easy to use graphical interface and access to several different analysis methods that allow the user to determine essentiality within a single condition as well as between two conditions.


            If you need to cite this tool, please use the following reference:

            DeJesus, M.A., Ambadipudi, C., Baker, R., Sassetti, C., and Ioerger, T.R. (2015). TRANSIT - a Software Tool for Himar1 TnSeq Analysis. PLOS Computational Biology, 11(10):e1004401


        """.replace("\n            ","\n")

        licence = """
            TRANSIT is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License.


            TRANSIT is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            GNU General Public License for more details.

            You should have received a copy of the GNU General Public License
            along with TRANSIT.  If not, see <http://www.gnu.org/licenses/>.
        """.replace("\n            ", "\n")

        info = wx.adv.AboutDialogInfo()
        info.SetIcon(images.transit_logo2.GetIcon())
        # images.transit_logo2.GetImage().ConvertToBitmap()
        info.SetName("TRANSIT")
        info.SetVersion(pytransit.__version__)
        info.SetDescription(description)
        info.SetCopyright("(C) 2015\n Michael A. DeJesus\nThomas R. Ioerger")
        info.SetWebSite("http://saclab.tamu.edu/essentiality/transit/")
        info.SetLicence(licence)
        info.AddDeveloper("Michael A. DeJesus")
        info.AddDeveloper("Thomas R. Ioerger")
        info.AddDeveloper("Chaitra Ambadipudi")
        info.AddDeveloper("Richard Baker")
        info.AddDeveloper("Christopher Sassetti")
        info.AddDeveloper("Eric Nelson")
        wx.adv.AboutBox(info)

    def documentationFunc(self, event):

        filepath = "http://saclab.tamu.edu/essentiality/transit/transit.html"
        output = ""
        error = ""
        try:
            if sys.platform.startswith("darwin"):
                OS = "OSX"
                # subprocess.call(('open', filepath))
                args = ["open", filepath]
                output, error = subprocess.Popen(
                    args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()
            elif os.name == "nt":
                OS = "Windows"
                os.startfile(filepath)
            elif os.name == "posix":
                OS = "Linux"
                args = ["xdg-open", filepath]
                output, error = subprocess.Popen(
                    args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()
                if "not found" in error:
                    args = ["exo-open", filepath]
                    output, error = subprocess.Popen(
                        args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    ).communicate()

        except Exception as e:
            error_text = """Error occurred opening documentation URL.\nYour browser or OS may not be configured correctly."""
            transit_tools.ShowError(MSG=error_text)
            traceback.print_exc()
    

    def MethodSelectFunc(self, selected_name, test=""):
        # If empty is selected
        if selected_name == "[Choose Method]":
            method_wrap_width = 250
            self.HideAllOptions()
            self.methodInfoText.SetLabel(u"Instructions")
            self.methodInstructions.Show()
            self.methodInstructions.SetLabel(self.instructions_text)
            self.methodInstructions.Wrap(method_wrap_width)
            self.methodShortText.Hide()
            self.methodLongText.Hide()
            self.methodTnText.Hide()
            self.methodDescText.Hide()

            self.method_choice = ""
        else:
            self.ShowGlobalOptions()
            self.methodSizerText.Show()

            # Get selected Method and hide Others
            for name in methods:
                methods[name].gui.Hide()
                methods[name].gui.GlobalHide()
                methods[name].gui.GlobalDisable()

                if methods[name].fullname() == selected_name:
                    matched_name = name

            if matched_name in methods:
                name = matched_name
                self.methodInfoText.SetLabel("%s" % methods[name].long_name)

                self.methodTnText.Show()
                self.methodTnText.SetLabel(methods[name].getTransposonsText())
                self.methodTnText.Wrap(250)

                self.methodDescText.Show()
                self.methodDescText.SetLabel(methods[name].getDescriptionText())
                self.methodDescText.Wrap(250)
                self.methodInstructions.SetLabel(" ")
                methods[name].gui.Show()
                methods[name].gui.Show()
                methods[name].gui.GlobalEnable()
                self.statusBar.SetStatusText("[%s]" % methods[name].short_name)

            self.ShowProgressSection()
            self.method_choice = selected_name

        self.Layout()
        if self.verbose:
            transit_tools.transit_message("Selected Method: %s" % (selected_name))

    def ExportSelectFunc(self, selected_name, test=""):
        # X = self.methodChoice.GetCurrentSelection()
        # selected_name = self.methodChoice.GetString(X)

        if self.verbose:
            transit_tools.transit_message(
                "Selected Export Method: %s" % (selected_name)
            )

        for name in export_methods:
            if export_methods[name].label == selected_name:
                methodobj = export_methods[name].method
                try:
                    M = methodobj.fromGUI(self)
                    if M:
                        thread = threading.Thread(target=M.Run())
                        thread.setDaemon(True)
                        thread.start()
                except Exception as e:
                    transit_tools.transit_message("Error: %s" % str(e))
                    traceback.print_exc()

    def ConvertSelectFunc(self, selected_name, test=""):
        annotationpath = self.annotation

        for name in convert_methods:
            if convert_methods[name].label == selected_name:
                methodobj = convert_methods[name].method
                try:
                    M = methodobj.fromGUI(self)
                    if M:
                        thread = threading.Thread(target=M.Run())
                        thread.setDaemon(True)
                        thread.start()
                except Exception as e:
                    transit_tools.transit_message("Error: %s" % str(e))
                    traceback.print_exc()

    def updateGraphChoices(self, dataset_type):

        empty_action = "[Choose Action]"
        if dataset_type == "Gumbel":
            choices = [empty_action, "Plot Ranked Probability of Essentiality"]
        elif dataset_type == "Binomial":
            choices = [empty_action, "Plot Ranked Probability of Essentiality"]
        elif dataset_type == "HMM - Sites":
            choices = [empty_action]
        elif dataset_type == "HMM - Genes":
            choices = [empty_action]
        elif dataset_type == "Resampling":
            choices = [
                empty_action,
                "Create a Volcano Plot",
                "Plot Histogram of logFC of Gene Counts",
            ]
        elif dataset_type == "DE-HMM - Sites":
            choices = [empty_action, "Recreate Sites File"]
        elif dataset_type == "DE-HMM - Segments":
            choices = [empty_action]
        else:
            choices = [empty_action]

        self.fileActionChoice.SetItems(choices)
        self.fileActionChoice.SetSelection(0)

    def graphGeneCounts(self, dataset_name, dataset_type, dataset_path):
        try:
            if dataset_type == "Resampling":
                X = []
                for line in open(dataset_path):
                    if line.startswith("#"):
                        continue
                    tmp = line.strip().split("\t")
                    try:
                        log2FC = float(tmp[-3])
                    except:
                        log2FC = 0
                    X.append(log2FC)

                n, bins, patches = plt.hist(
                    X, density=1, facecolor="c", alpha=0.75, bins=100
                )
                plt.xlabel("log2 FC - Total Gene Counts")
                plt.ylabel("Probability")
                plt.title(
                    "Histogram of log2 Fold Change for Total Normalized Counts within Genes"
                )
                plt.axvline(0, color="r", linestyle="dashed", linewidth=3)
                plt.grid(True)
                plt.show()
            else:
                transit_tools.ShowError(
                    MSG="Need to select a 'Resampling' results file for this type of plot."
                )

        except Exception as e:
            transit_tools.transit_message("Error occurred creating plot: %s" % str(e))
            traceback.print_exc()

    def graphVolcanoPlot(self, dataset_name, dataset_type, dataset_path):
        try:
            if dataset_type == "Resampling":
                X = []
                Y = []
                header = []
                qval_list = []
                bad = []
                col_logFC = -6
                col_pval = -2
                col_qval = -1
                ii = 0
                for line in open(dataset_path):
                    if line.startswith("#"):
                        tmp = line.split("\t")
                        temp_col_logfc = [
                            i
                            for (i, x) in enumerate(tmp)
                            if "logfc" in x.lower()
                            or "log-fc" in x.lower()
                            or "log2fc" in x.lower()
                        ]
                        temp_col_pval = [
                            i
                            for (i, x) in enumerate(tmp)
                            if ("pval" in x.lower() or "p-val" in x.lower())
                            and "adj" not in x.lower()
                        ]
                        if temp_col_logfc:
                            col_logFC = temp_col_logfc[-1]
                        if temp_col_pval:
                            col_pval = temp_col_pval[-1]
                        continue

                    tmp = line.strip().split("\t")
                    try:
                        log10qval = -math.log(float(tmp[col_pval].strip()), 10)
                    except ValueError as e:
                        bad.append(ii)
                        log10qval = 0

                    log2FC = float(tmp[col_logFC])

                    qval_list.append(
                        (float(tmp[col_qval]), float(tmp[col_pval].strip()))
                    )
                    X.append(log2FC)
                    Y.append(log10qval)
                    ii += 1
                count = 0
                threshold = 0.00001
                backup_thresh = 0.00001
                qval_list.sort()
                for (q, p) in qval_list:
                    backup_thresh = p
                    if q > 0.05:
                        break
                    threshold = p
                    count += 1

                if threshold == 0:
                    threshold = backup_thresh
                for ii in bad:
                    Y[ii] = max(Y)
                plt.plot(X, Y, "bo")
                plt.axhline(
                    -math.log(threshold, 10), color="r", linestyle="dashed", linewidth=3
                )
                plt.xlabel("Log Fold Change (base 2)")
                plt.ylabel("-Log p-value (base 10)")
                plt.suptitle("Resampling - Volcano plot")
                plt.title("Adjusted threshold (red line): %1.8f" % threshold)
                plt.show()
            else:
                transit_tools.ShowError(
                    MSG="Need to select a 'Resampling' results file for this type of plot."
                )

        except Exception as e:
            print("Error occurred creating plot:", str(e))

    def graphRankedZbar(self, dataset_name, dataset_type, dataset_path):
        try:
            X = []
            Y = []
            for line in open(dataset_path):
                if line.startswith("#"):
                    continue
                tmp = line.strip().split("\t")
                try:
                    # log2FC = math.log(float(tmp[6])/float(tmp[5]),2)
                    zbar = float(tmp[-2])
                except:
                    zbar = 0
                if zbar >= 0:
                    Y.append(zbar)

            Y.sort()
            index = range(1, len(Y) + 1)
            plt.plot(index, Y, "bo")
            plt.xlabel("Rank of Gene According to Probability of Essentiality")
            plt.ylabel("Probability of Essentiality")
            plt.title("Ranked Probability of Essentiality")
            plt.show()

        except Exception as e:
            print("Error occurred creating plot:", str(e))

    def LoessPrevFunc(self, event):
        datasets_selected = self.ctrlSelected() + self.expSelected()
        if not datasets_selected:
            transit_tools.ShowError(
                MSG="Need to select at least one control or experimental dataset."
            )
            return

        data, position = tnseq_tools.get_data(datasets_selected)
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

    def choseMethodsMenu(self, selected_name, event):
        if self.verbose:
            transit_tools.transit_message("Selected Method: %s" % (selected_name))
        self.MethodSelectFunc(selected_name)

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
            transit_tools.transit_message(
                "Selected the '%s' normalization method" % dlg.GetStringSelection()
            )

        dlg.Destroy()
        return dlg.GetStringSelection()

    def annotationPT_to_GFF3(self, event):
        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".gff3"
        # defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()
        outputPath = self.SaveFile(defaultDir, defaultFile)

        ORGANISM = transit_tools.fetch_name(annotationpath)
        if not annotationpath:
            transit_tools.ShowError("Error: No annotation file selected.")

        elif outputPath:
            if self.verbose:
                transit_tools.transit_message(
                    "Converting annotation file from prot_table format to GFF3 format"
                )
            year = time.localtime().tm_year
            month = time.localtime().tm_mon
            day = time.localtime().tm_mday

            output = open(outputPath, "w")
            output.write("##gff-version 3\n")
            output.write("##converted to IGV with TRANSIT.\n")
            output.write("##date %d-%d-%d\n" % (year, month, day))
            output.write("##Type DNA %s\n" % ORGANISM)

            for line in open(annotationpath):
                if line.startswith("#"):
                    continue
                tmp = line.strip().split("\t")
                desc = tmp[0]
                start = int(tmp[1])
                end = int(tmp[2])
                strand = tmp[3]
                length = tmp[4]
                name = tmp[7]
                orf = tmp[8]
                ID = name
                desc.replace("%", "%25").replace(";", "%3B").replace(
                    "=", "%3D"
                ).replace(",", "%2C")
                output.write(
                    "%s\tRefSeq\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s;Alias=%s;locus_tag=%s;desc=%s\n"
                    % (ORGANISM, start, end, strand, orf, ID, orf, orf, desc)
                )

            output.close()
            if self.verbose:
                transit_tools.transit_message("Finished conversion")

    def annotationPT_to_PTT(self, event):

        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".ptt.table"
        # defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        datasets = self.ctrlSelected() + self.expSelected()
        if not annotationpath:
            transit_tools.ShowError("Error: No annotation file selected.")
        elif not datasets:
            transit_tools.ShowError(
                "Error: Please add a .wig dataset, to determine TA sites."
            )
        else:

            outputPath = self.SaveFile(defaultDir, defaultFile)
            if not outputPath:
                return
            if self.verbose:
                transit_tools.transit_message(
                    "Converting annotation file from prot_table format to PTT format"
                )
            (data, position) = tnseq_tools.get_data(datasets)
            orf2info = transit_tools.get_gene_info(annotationpath)
            hash = transit_tools.get_pos_hash(annotationpath)
            (orf2reads, orf2pos) = tnseq_tools.get_gene_reads(
                hash, data, position, orf2info
            )

            output = open(outputPath, "w")
            output.write("geneID\tstart\tend\tstrand\tTA coordinates\n")
            for line in open(annotationpath):
                if line.startswith("#"):
                    continue
                tmp = line.strip().split("\t")
                orf = tmp[8]
                name = tmp[7]
                desc = tmp[0]
                start = int(tmp[1])
                end = int(tmp[2])
                strand = tmp[3]
                ta_str = "no TAs"
                if orf in orf2pos:
                    ta_str = "\t".join([str(int(ta)) for ta in orf2pos[orf]])
                output.write("%s\t%s\t%s\t%s\t%s\n" % (orf, start, end, strand, ta_str))
            output.close()
            if self.verbose:
                transit_tools.transit_message("Finished conversion")

    def annotationPTT_to_PT(self, event):

        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".prot_table"
        # defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        datasets = self.ctrlSelected() + self.expSelected()
        if not annotationpath:
            transit_tools.ShowError("Error: No annotation file selected.")
        # elif not datasets:
        #    transit_tools.ShowError("Error: Please add a .wig dataset, to determine TA sites.")
        else:

            outputPath = self.SaveFile(defaultDir, defaultFile)
            if not outputPath:
                return
            if self.verbose:
                transit_tools.transit_message(
                    "Converting annotation file from PTT format to prot_table format"
                )
            # (data, position) = tnseq_tools.get_data(datasets)
            # orf2info = transit_tools.get_gene_info(annotationpath)
            # hash = transit_tools.get_pos_hash(annotationpath)
            # (orf2reads, orf2pos) = tnseq_tools.get_gene_reads(hash, data, position, orf2info)

            output = open(outputPath, "w")
            # output.write("geneID\tstart\tend\tstrand\tTA coordinates\n")
            for line in open(annotationpath):
                if line.startswith("#"):
                    continue
                if line.startswith("geneID"):
                    continue
                tmp = line.strip().split("\t")
                orf = tmp[0]
                if orf == "intergenic":
                    continue
                name = "-"
                desc = "-"
                start = int(tmp[1])
                end = int(tmp[2])
                length = ((end - start + 1) / 3) - 1
                strand = tmp[3]
                someID = "-"
                someID2 = "-"
                COG = "-"
                output.write(
                    "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n"
                    % (
                        desc,
                        start,
                        end,
                        strand,
                        length,
                        someID,
                        someID2,
                        name,
                        orf,
                        COG,
                    )
                )
            output.close()
            if self.verbose:
                transit_tools.transit_message("Finished conversion")

    def annotationGFF3_to_PT(self, event):

        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".prot_table"
        # defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        datasets = self.ctrlSelected() + self.expSelected()
        if not annotationpath:
            transit_tools.ShowError("Error: No annotation file selected.")
        else:
            outputPath = self.SaveFile(defaultDir, defaultFile)
            if not outputPath:
                return
            if self.verbose:
                transit_tools.transit_message(
                    "Converting annotation file from GFF3 format to prot_table format"
                )

            output = open(outputPath, "w")
            for line in open(annotationpath):
                if line.startswith("#"):
                    continue
                tmp = line.strip().split("\t")
                chr = tmp[0]
                type = tmp[2]
                start = int(tmp[3])
                end = int(tmp[4])
                length = ((end - start + 1) / 3) - 1
                strand = tmp[6]
                features = dict([tuple(f.split("=")) for f in tmp[8].split(";")])
                if "ID" not in features:
                    continue
                orf = features["ID"]
                name = features.get("Name", "-")
                if name == "-":
                    name = features.get("name", "-")

                desc = features.get("Description", "-")
                if desc == "-":
                    desc = features.get("description", "-")
                if desc == "-":
                    desc = features.get("Desc", "-")
                if desc == "-":
                    desc = features.get("desc", "-")

                someID = "-"
                someID2 = "-"
                COG = "-"
                output.write(
                    "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n"
                    % (
                        desc,
                        start,
                        end,
                        strand,
                        length,
                        someID,
                        someID2,
                        name,
                        orf,
                        COG,
                    )
                )
            output.close()
            if self.verbose:
                transit_tools.transit_message("Finished conversion")

    def RunMethod(self, event):
        # FLORF
        # X = self.methodChoice.GetCurrentSelection()
        # selected_name = self.methodChoice.GetString(X)
        selected_name = self.method_choice
        for name in methods:
            if methods[name].fullname() == selected_name:
                methodobj = methods[name].method
        try:
            M = methodobj.fromGUI(self)
            if M:
                thread = threading.Thread(target=M.Run())
                thread.setDaemon(True)
                thread.start()
        except Exception as e:
            transit_tools.transit_message("Error: %s" % str(e))
            traceback.print_exc()
