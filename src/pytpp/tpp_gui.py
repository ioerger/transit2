#!/usr/bin/env python

# Copyright 2017.
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

import sys
import glob
import os
import time
import math
import re
import shutil
import platform
import gzip

try:
    import wx
    import wx.lib.filebrowsebutton

    HAS_WX = True
except Exception as e:
    HAS_WX = False

from pytpp.tpp_tools import *

if HAS_WX:

    class TPPIcon(wx.StaticBitmap):
        def __init__(self, panel, flag, bit_map, tooltip=""):
            wx.StaticBitmap.__init__(self, panel, flag, bit_map)
            tp = wx.ToolTip(tooltip)
            self.SetToolTip(tp)

    class MyForm(wx.Frame):
        def __init__(self, vars):
            self.vars = vars
            initialize_globals(self.vars)

            wx.Frame.__init__(
                self, None, wx.ID_ANY, "TPP: Tn-Seq PreProcessor"
            )  # v%s" % vars.version
            # Add a panel so it looks the correct on all platforms

            panel = wx.ScrolledWindow(
                self,
                wx.ID_ANY,
                wx.DefaultPosition,
                wx.Size(-1, -1),
                wx.HSCROLL | wx.VSCROLL,
            )
            panel.SetScrollRate(5, 5)
            panel.SetMaxSize(wx.Size(-1, 1000))

            sizer = wx.BoxSizer(wx.VERTICAL)

            self.listCtrl = None
            self.InitMenu()
            self.InitFiles(panel, sizer)

            buttonrow = wx.BoxSizer(wx.HORIZONTAL)

            btn = wx.Button(panel, label="Start")
            btn.Bind(wx.EVT_BUTTON, self.map_reads)
            buttonrow.Add(btn, 0, 0, 0, 10)

            btn = wx.Button(panel, label="Quit")
            btn.Bind(wx.EVT_BUTTON, self.OnQuit)
            buttonrow.Add(btn, 0, 0, 0, 10)
            sizer.Add(buttonrow, 0, 0, 0)

            self.InitList(panel, sizer)

            panel.SetSizer(sizer)
            # self.SetSize((1305, 700))
            self.SetSize((900, 750))
            # self.SetTitle('Simple menu')
            self.Centre()
            # self.Show(True)

            self.pid = None

        #

        def InitFiles(self, panel, sizer):
            vars = self.vars
            # Define
            bit_map = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (16, 16))

            # REFERENCE
            sizer3 = wx.BoxSizer(wx.HORIZONTAL)
            label3 = wx.StaticText(
                panel,
                label="Choose a reference genome (FASTA) (REQUIRED):",
                size=(330, -1),
            )
            sizer3.Add(label3, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.picker3 = wx.lib.filebrowsebutton.FileBrowseButton(
                panel,
                id=wx.ID_ANY,
                dialogTitle="Please select the reference genome",
                fileMode=wx.FD_OPEN,
                fileMask="*.fna;*.fasta;*.fa",
                size=(400, 30),
                startDirectory=os.path.dirname(vars.ref),
                initialValue=vars.ref,
                labelText="",
            )
            sizer3.Add(self.picker3, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            sizer3.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Select a reference genome in FASTA format (can be a multi-contig fasta file).",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer3.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer3, 0, wx.EXPAND, 0)

            # REPLICON ID NAMES
            sizer_replicon_ids = wx.BoxSizer(wx.HORIZONTAL)
            label_replicon_ids = wx.StaticText(
                panel,
                label="ID names for each replicon: \n(if genome has multiple contigs)",
                size=(340, -1),
            )
            sizer_replicon_ids.Add(label_replicon_ids, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.replicon_ids = wx.TextCtrl(
                panel, value=vars.replicon_ids, size=(400, 30)
            )
            sizer_replicon_ids.Add(
                self.replicon_ids, proportion=1.0, flag=wx.EXPAND | wx.ALL, border=5
            )
            sizer_replicon_ids.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Specify names of each contig within the reference genome separated by commas (if using wig_gb_to_csv.py you must use the contig names in the Genbank file).  Only required if there are multiple contigs; can leave blank if there is just one sequence.\nEnter 'auto' for autogenerated ids.",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer_replicon_ids.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer_replicon_ids, 0, wx.EXPAND, 0)

            # READS 1
            sizer1 = wx.BoxSizer(wx.HORIZONTAL)
            label1 = wx.StaticText(
                panel,
                label="Choose the Fastq file for read 1 (REQUIRED):",
                size=(330, -1),
            )
            sizer1.Add(label1, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.picker1 = wx.lib.filebrowsebutton.FileBrowseButton(
                panel,
                id=wx.ID_ANY,
                dialogTitle="Please select the .fastq file for read 1",
                fileMode=wx.FD_OPEN,
                fileMask="*.fastq;*.fq;*.reads;*.fasta;*.fa;*.fastq.gz",
                size=(400, 30),
                startDirectory=os.path.dirname(vars.fq1),
                initialValue=vars.fq1,
                labelText="",
                changeCallback=self.OnChanged2,
            )
            sizer1.Add(self.picker1, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            sizer1.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Select a file containing the reads in .FASTQ (or compressed FASTQ) format.",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer1.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer1, 0, wx.EXPAND, 0)

            # READS 2
            sizer2 = wx.BoxSizer(wx.HORIZONTAL)
            label2 = wx.StaticText(
                panel, label="Choose the Fastq file for read 2:", size=(330, -1)
            )
            sizer2.Add(label2, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.picker2 = wx.lib.filebrowsebutton.FileBrowseButton(
                panel,
                id=wx.ID_ANY,
                dialogTitle="Please select the .fastq file for read 2",
                fileMode=wx.FD_OPEN,
                fileMask="*.fastq;*.fq;*.reads;*.fasta;*.fa;*.fastq.gz",
                size=(400, 30),
                startDirectory=os.path.dirname(vars.fq2),
                initialValue=vars.fq2,
                labelText="",
                changeCallback=self.OnChanged2,
            )
            sizer2.Add(self.picker2, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            sizer2.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Select a file containing the pair-end reads in .FASTQ (or compressed FASTQ) format. Optional.",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer2.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer2, 0, wx.EXPAND, 0)

            # OUTPUT PREFIX
            sizer5 = wx.BoxSizer(wx.HORIZONTAL)
            label5 = wx.StaticText(
                panel,
                label="Prefix to use for output filenames (REQUIRED):",
                size=(340, -1),
            )
            sizer5.Add(label5, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.base = wx.TextCtrl(panel, value=vars.base, size=(400, 30))
            sizer5.Add(self.base, proportion=1.0, flag=wx.EXPAND | wx.ALL, border=5)
            sizer5.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Select a prefix that will be used when writing output files",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer5.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer5, 0, wx.EXPAND, 0)

            # PROTOCOL
            sizer_protocol = wx.BoxSizer(wx.HORIZONTAL)
            label_protocol = wx.StaticText(
                panel, label="Protocol used:", size=(340, -1)
            )
            sizer_protocol.Add(label_protocol, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.protocol = wx.ComboBox(
                panel, choices=["Sassetti", "Mme1", "Tn5"], size=(400, 30)
            )
            self.protocol.SetStringSelection(vars.protocol)
            sizer_protocol.Add(
                self.protocol, proportion=1, flag=wx.EXPAND | wx.ALL, border=5
            )
            protocol_tooltip_text = """Select which protocol used to prepare the sequencing samples. Default values will populate the other fields.

The Sassetti protocol generally assumes the reads include the primer prefix and part of the transposon sequence, followed by genomic sequence. It also assumes reads are sequenced in the forward direction.  Barcodes are in read 2, along with genomic DNA from the other end of the fragment.

The Mme1 protocol generally assumes reads do NOT include the primer prefix, and that the reads are sequenced in the reverse direction"""
            sizer_protocol.Add(
                TPPIcon(panel, wx.ID_ANY, bit_map, protocol_tooltip_text),
                flag=wx.CENTER,
                border=0,
            )
            sizer_protocol.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer_protocol, 0, wx.EXPAND, 0)

            self.Bind(
                wx.EVT_COMBOBOX, self.OnProtocolSelection, id=self.protocol.GetId()
            )

            # TRANSPOSON
            sizer8 = wx.BoxSizer(wx.HORIZONTAL)
            label8 = wx.StaticText(panel, label="Transposon used:", size=(340, -1))
            sizer8.Add(label8, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.transposon = wx.ComboBox(
                panel,
                choices=["Himar1", "Tn5", "pre-trimmed", "[Custom]"],
                size=(400, 30),
            )
            self.transposon.SetStringSelection(vars.transposon)
            sizer8.Add(self.transposon, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            sizer8.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Select the transposon used to construct the TnSeq libraries. This will automatically populate the primer prefix field. Select custom to specify your own sequence.",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer8.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer8, 0, wx.EXPAND, 0)

            # PRIMER SEQUENCE
            sizer4 = wx.BoxSizer(wx.HORIZONTAL)
            label4 = wx.StaticText(panel, label="Primer sequence:", size=(340, -1))
            sizer4.Add(label4, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.prefix = wx.TextCtrl(panel, value=str(vars.prefix), size=(400, 30))
            sizer4.Add(self.prefix, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            sizer4.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "If present in the reads, specify the primer sequence. If it has been stripped away already, leave this field empty.",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer4.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer4, 0, wx.EXPAND, 0)

            self.Bind(
                wx.EVT_COMBOBOX, self.OnTransposonSelection, id=self.transposon.GetId()
            )
            self.prefix.Bind(wx.EVT_TEXT, self.OnChangePrimerPrefix)

            # MAX READS
            sizer6 = wx.BoxSizer(wx.HORIZONTAL)
            label6 = wx.StaticText(
                panel, label="Max reads (leave blank to use all):", size=(340, -1)
            )
            sizer6.Add(label6, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.maxreads = wx.TextCtrl(
                panel, value=str(vars.maxreads), size=(150, 30)
            )  # or "" if not defined? can't write to tpp.cfg
            sizer6.Add(self.maxreads, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            sizer6.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Maximum reads to use from the reads files. Useful for running only a portion of very large number of reads. Leave blank to use all the reads.",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer6.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer6, 0, wx.EXPAND, 0)

            # MISMATCHES
            sizer7 = wx.BoxSizer(wx.HORIZONTAL)
            label7 = wx.StaticText(
                panel, label="Mismatches allowed in Tn prefix:", size=(340, -1)
            )
            sizer7.Add(label7, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.mismatches = wx.TextCtrl(panel, value=str(vars.mm1), size=(150, 30))
            sizer7.Add(self.mismatches, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            sizer7.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Number of mismatches allowed in the tn-prefix before discarding the read.",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer7.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer7, 0, wx.EXPAND, 0)

            # PRIMER_START_WINDOW
            sizer_primer_start = wx.BoxSizer(wx.HORIZONTAL)
            label_primer_start = wx.StaticText(
                panel,
                label="Start of window to look for prefix (Tn terminus):",
                size=(340, -1),
            )
            sizer_primer_start.Add(label_primer_start, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            primer_start_window = "%s,%s" % (
                vars.primer_start_window[0],
                vars.primer_start_window[1],
            )
            self.primer_start = wx.TextCtrl(
                panel, value=primer_start_window, size=(150, 30)
            )
            sizer_primer_start.Add(
                self.primer_start, proportion=1, flag=wx.EXPAND | wx.ALL, border=5
            )
            sizer_primer_start.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Region in read 1 to search for start of prefix seq (i.e. end of transposon).",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer_primer_start.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer_primer_start, 0, wx.EXPAND, 0)

            #            # WINDOW SIZE                                 # [RJ] This block is to add the acceptance of a set window size for setting P,Q parameters
            #            sizer_window_size = wx.BoxSizer(wx.HORIZONTAL)
            #            label_window_size = wx.StaticText(panel, label='Window size for Tn prefix in read:', size=(340,-1))
            #            sizer_window_size.Add(label_window_size,0,wx.ALIGN_CENTER_VERTICAL,0)
            #            self.window_size = wx.TextCtrl(panel,value=str(vars.window_size),size=(150,30))
            #            sizer_window_size.Add(self.window_size, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
            #            sizer_window_size.Add(TPPIcon(panel, wx.ID_ANY, bit_map, "Window size for extract_staggered() to look for start of Tn prefix."), flag=wx.CENTER, border=0)
            #            sizer_window_size.Add((10, 1), 0, wx.EXPAND)
            #            sizer.Add(sizer_window_size,0,wx.EXPAND,0)

            # BWA
            sizer0 = wx.BoxSizer(wx.HORIZONTAL)
            label0 = wx.StaticText(
                panel, label="BWA executable (REQUIRED):", size=(330, -1)
            )
            sizer0.Add(label0, 0, wx.ALIGN_CENTER_VERTICAL, 0)

            self.picker0 = wx.lib.filebrowsebutton.FileBrowseButton(
                panel,
                id=wx.ID_ANY,
                size=(400, 30),
                dialogTitle="Path to BWA",
                fileMode=wx.FD_OPEN,
                fileMask="bwa*",
                startDirectory=os.path.dirname(vars.bwa),
                initialValue=vars.bwa,
                labelText="",
            )
            sizer0.Add(self.picker0, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            sizer0.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Specify a path to the BWA executable (including the executable).",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer0.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer0, 0, wx.EXPAND, 0)

            self.bwa_alg = wx.ComboBox(
                panel,
                choices=["use algorithm 'aln'", "use algorithm 'mem'"],
                size=(200, 30),
            )
            if vars.bwa_alg == "aln":
                self.bwa_alg.SetSelection(0)
            else:
                self.bwa_alg.SetSelection(1)  # default
            sizer0.Add(
                self.bwa_alg, proportion=0.5, flag=wx.EXPAND | wx.ALL, border=5
            )  ##
            self.bwa_alg.Bind(
                wx.EVT_COMBOBOX, self.OnBwaAlgSelection, id=self.bwa_alg.GetId()
            )
            sizer0.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "'mem' is considered to do a better job at mapping reads, but 'aln' is available as an alternative.",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer0.Add((10, 1), 0, wx.EXPAND)
            # sizer.Add(sizer0,0,wx.EXPAND,0)

            # BWA FLAGS
            sizer8 = wx.BoxSizer(wx.HORIZONTAL)
            label8 = wx.StaticText(panel, label="BWA flags (Optional)", size=(340, -1))
            sizer8.Add(label8, 0, wx.ALIGN_CENTER_VERTICAL, 0)
            self.flags = wx.TextCtrl(panel, value=vars.flags, size=(400, 30))
            sizer8.Add(self.flags, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            sizer8.Add(
                TPPIcon(
                    panel,
                    wx.ID_ANY,
                    bit_map,
                    "Use this textbox to enter any desired flags for the BWA alignment. For example, to limit the number of mismatches to 1, type: -k 1. See the BWA documentation for all possible flags.",
                ),
                flag=wx.CENTER,
                border=0,
            )
            sizer8.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer8, 0, wx.EXPAND, 0)

            # BARSEQ CATALOG
            sizer9 = wx.BoxSizer(wx.HORIZONTAL)
            label9 = wx.StaticText(panel, label="BarSeq Catalog file:", size=(120, -1))
            sizer9.Add(label9, 0, wx.ALIGN_CENTER_VERTICAL, 0)

            self.barseq_select = wx.ComboBox(
                panel,
                choices=[
                    "this is not a Barseq dataset",
                    "read catalog file",
                    "write catalog file",
                ],
                size=(200, 30),
            )  ## # does a BoxSizer use wx.HORIZONTAL, not wx.EXPAND?
            self.barseq_select.SetSelection(0)
            sizer9.Add(
                self.barseq_select, proportion=0.5, flag=wx.EXPAND | wx.ALL, border=5
            )  ##

            self.picker9 = wx.lib.filebrowsebutton.FileBrowseButton(
                panel,
                id=wx.ID_ANY,
                dialogTitle="Please select the Barseq catalog filename",
                fileMode=wx.FD_OPEN,
                size=(400, 30),
                startDirectory=os.path.dirname(vars.fq2),
                initialValue="",
                labelText="",
            )  # no need for this: changeCallback=self.OnChanged9 ; initialValue set below ; no file mask
            sizer9.Add(self.picker9, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)

            if vars.barseq_catalog_in != None:
                self.barseq_select.SetSelection(1)
                self.picker9.SetValue(vars.barseq_catalog_in)
            if vars.barseq_catalog_out != None:
                self.barseq_select.SetSelection(2)
                self.picker9.SetValue(vars.barseq_catalog_out)

            sizer9.Add(
                TPPIcon(panel, wx.ID_ANY, bit_map, "Select a filename for BarSeq catalog."),
                flag=wx.CENTER,
                border=0,
            )
            sizer9.Add((10, 1), 0, wx.EXPAND)
            sizer.Add(sizer9, 0, wx.EXPAND, 0)

        #

        def OnBwaAlgSelection(self, event):
            if "aln" in self.bwa_alg.GetValue():
                self.vars.bwa_alg = "aln"
            elif "mem" in self.bwa_alg.GetValue():
                self.vars.bwa_alg = "mem"
            else:
                self.vars.bwa_alg = "[Custom]"

        #

        def OnTransposonSelection(self, event):
            if self.transposon.GetValue() == "Tn5":
                self.prefix.SetValue("TAAGAGACAG")
                self.transposon.SetStringSelection("Tn5")
                self.vars.transposon = "Tn5"
            elif self.transposon.GetValue() == "Himar1":
                self.prefix.SetValue("ACTTATCAGCCAACCTGTTA")
                self.transposon.SetStringSelection("Himar1")
                self.vars.transposon = "Himar1"
            elif self.transposon.GetValue() == "pre-trimmed":
                self.transposon.SetValue("pre-trimmed")
                self.transposon.SetStringSelection("pre-trimmed")
                self.vars.transposon = "pre-trimmed"
                self.prefix.SetValue('""')
            else:
                self.transposon.SetValue("[Custom]")
                self.transposon.SetStringSelection("[Custom]")
                self.vars.transposon = "[Custom]"

        #

        def OnProtocolSelection(self, event):
            self.vars.transposon = self.protocol.GetValue()
            if self.protocol.GetValue() == "Tn5":
                self.prefix.SetValue("TAAGAGACAG")
                self.transposon.SetStringSelection("Tn5")
                self.vars.transposon = "Tn5"
            elif self.protocol.GetValue() == "Sassetti":
                self.prefix.SetValue("ACTTATCAGCCAACCTGTTA")
                self.transposon.SetStringSelection("Himar1")
                self.vars.transposon = "Himar1"
            elif self.protocol.GetValue() == "Mme1":
                self.prefix.SetValue('""')
                self.transposon.SetStringSelection("pre-trimmed")
                self.vars.transposon = ""

        #

        def OnChanged(self, str_path):
            print("changed")
            value = os.path.basename(str_path).split(".")[0]
            if "_R1" in value or "_R2":
                value = value.split("_")[0]
            # self.base.SetValue(value)

        #

        def OnChanged2(self, event):
            value2 = os.path.basename(self.picker2.GetValue()).split(".")[0]
            value1 = os.path.basename(self.picker1.GetValue()).split(".")[0]
            value = os.path.commonprefix([value1, value2])
            # self.base.SetValue(value)
            # self.base.Refresh()

        #

        def OnChangePrimerPrefix(self, event):
            # self.transposon.SetValue("[Custom]")
            pass

        #

        def InitList(self, panel, sizer):
            self.listCtrl = wx.ListCtrl(
                panel, size=(500, 210), style=wx.LC_REPORT | wx.BORDER_SUNKEN
            )
            self.listCtrl.InsertColumn(0, "Dataset (*.tn_stats)", width=300)
            self.listCtrl.InsertColumn(
                1, "total reads", wx.LIST_FORMAT_RIGHT, width=125
            )
            self.listCtrl.InsertColumn(2, "Tn prefix", wx.LIST_FORMAT_RIGHT, width=125)
            self.listCtrl.InsertColumn(3, "R1_mapped", wx.LIST_FORMAT_RIGHT, width=90)
            self.listCtrl.InsertColumn(4, "R2_mapped", wx.LIST_FORMAT_RIGHT, width=90)
            self.listCtrl.InsertColumn(
                5, "mapped\nreads", wx.LIST_FORMAT_RIGHT, width=90
            )
            self.listCtrl.InsertColumn(
                6, "template\ncount", wx.LIST_FORMAT_RIGHT, width=90
            )
            self.listCtrl.InsertColumn(7, "TAs hit", wx.LIST_FORMAT_RIGHT, width=90)
            self.listCtrl.InsertColumn(
                8, "insertion\ndensity", wx.LIST_FORMAT_RIGHT, width=90
            )
            self.listCtrl.InsertColumn(9, "NZmean", wx.LIST_FORMAT_RIGHT, width=90)
            self.listCtrl.InsertColumn(10, "maxcount", wx.LIST_FORMAT_RIGHT, width=90)
            self.listCtrl.InsertColumn(11, "primer", wx.LIST_FORMAT_RIGHT, width=90)
            self.listCtrl.InsertColumn(12, "vector", wx.LIST_FORMAT_RIGHT, width=90)

            sizer.Add(self.listCtrl, 0, wx.ALL | wx.EXPAND, 10)

        #

        def InitMenu(self):
            menubar = wx.MenuBar()
            fileMenu = wx.Menu()

            quit_menuitem = fileMenu.Append(wx.ID_EXIT, "Quit", "Quit application")
            self.Bind(wx.EVT_MENU, self.OnQuit, quit_menuitem)

            menubar.Append(fileMenu, "&File")
            self.SetMenuBar(menubar)

        #

        def addNewDataset(self, event):
            dlg = wx.FileDialog(
                self,
                message="Choose a file",
                defaultDir=".",
                defaultFile="",
                wildcard="*.wig",
                style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
            )
            if dlg.ShowModal() == wx.ID_OK:
                paths = dlg.GetPaths()
                for path in paths:
                    print("analyzing dataset:", path)
                    analyze_dataset(path)
            dlg.Destroy()
            self.update_dataset_list()

        #

        def update_dataset_list(self):
            if self.listCtrl == None:
                return
            self.listCtrl.DeleteAllItems()
            self.index = 0
            datasets = []
            for fname in glob.glob("*.tn_stats"):
                filedate = os.path.getmtime(fname)
                datasets.append((filedate, fname))
            datasets.sort(reverse=True)
            for (filedate, fname) in datasets:
                stats = self.read_stats_file(fname)
                ntrim = stats.get("TGTTA_reads", "?")
                if ntrim == "?":
                    ntrim = stats.get("trimmed_reads", "?")
                vals = [
                    stats.get("total_reads", "?"),
                    ntrim,
                    stats.get("reads1_mapped", "?"),
                    stats.get("reads2_mapped", "?"),
                    stats.get("mapped_reads", "?"),
                    stats.get("template_count", "?"),
                    stats.get("TAs_hit", "?"),
                    stats.get("density", "?"),
                    stats.get("NZ_mean", "?"),
                    stats.get("max_count", "?"),
                    stats.get("primer_matches:", "?"),
                    stats.get("vector_matches:", "?"),
                ]

                dsname = "[%s] %s" % (
                    time.strftime("%m/%d/%y", time.localtime(filedate)),
                    fname[: fname.rfind(".")],
                )
                self.add_data(dsname, vals)

        #

        def read_stats_file(self, fname):
            stats = {}
            with open(fname) as file:
                for line in file:
                    w = line.rstrip().split()
                    val = ""
                    if len(w) > 2:
                        val = w[2]
                    stats[w[1]] = val
            return stats

        #

        def add_data(self, dataset, vals):
            self.listCtrl.InsertItem(self.index, dataset)
            for i in range(1, len(vals) + 1):
                self.listCtrl.SetItem(self.index, i, vals[i - 1])
            self.index += 1

        #

        def OnQuit(self, e):
            print("Quitting TPP.  Good bye.")
            self.vars.action = "quit"
            self.Close()
            return 0

        #

        def map_reads(self, event):

            # add bwa path, prefix
            bwapath = self.picker0.GetValue()
            fq1, fq2, ref, base, prefix, maxreads = (
                self.picker1.GetValue(),
                self.picker2.GetValue(),
                self.picker3.GetValue(),
                self.base.GetValue(),
                self.prefix.GetValue(),
                self.maxreads.GetValue(),
            )

            mm1 = self.mismatches.GetValue()
            try:
                mm1 = int(mm1)
            except Exception:
                mm1 = 1

            self.vars.flags = self.flags.GetValue()

            self.vars.transposon = self.transposon.GetStringSelection()
            self.vars.protocol = self.protocol.GetValue()

            self.vars.bwa = bwapath
            self.vars.fq1 = fq1
            self.vars.fq2 = fq2
            self.vars.ref = ref
            self.vars.base = base
            self.vars.mm1 = mm1
            self.vars.prefix = prefix

            # self.vars.window_size = int(self.window_size.GetValue())
            if "aln" in self.bwa_alg.GetValue():
                self.vars.bwa_alg = "aln"
            elif "mem" in self.bwa_alg.GetValue():
                self.vars.bwa_alg = "mem"
            self.vars.replicon_ids = self.replicon_ids.GetValue().split(",")

            v = self.primer_start.GetValue()
            if v != "":
                v = v.split(",")
                self.vars.primer_start_window = (int(v[0]), int(v[1]))

            if maxreads == "":
                self.vars.maxreads = -1
            else:
                self.vars.maxreads = int(maxreads)

            barseq_select = self.barseq_select.GetSelection()
            self.vars.barseq_catalog_in = self.vars.barseq_catalog_out = None
            if barseq_select == 1:
                self.vars.barseq_catalog_in = self.picker9.GetValue()
            if barseq_select == 2:
                self.vars.barseq_catalog_out = self.picker9.GetValue()

            self.vars.action = "start"
            self.Close()
            return 0
