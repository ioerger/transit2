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
from pytransit.core_data import SessionData, universal
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
        self.Bind(wx.EVT_TIMER, self.clearStatus, self.timer)

        
        self.SetIcon(images.transit_icon.GetIcon())

        self.workdir = os.getcwd()
        self.annotation = ""
        self.transposons = ["himar1", "tn5"]
        self.verbose = True
        
        self.statusBar = self.CreateStatusBar(1, wx.STB_SIZEGRIP, wx.ID_ANY)
        self.statusBar.SetStatusText("Welcome to TRANSIT")
        
        pub.subscribe(self.saveHistogram, "histogram")
        create_menu(self)
        
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
            transit_tools.log("Error Modifying Ctrl Lib String: %s" % e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def ctrlRemoveFunc(self, event):
        next = self.wig_table.GetNextSelected(-1)
        while next != -1:
            if self.verbose:
                transit_tools.log(
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
                transit_tools.log(
                    "Visualizing counts for: %s"
                    % ", ".join([transit_tools.fetch_name(d) for d in datasets])
                )
            viewWindow = trash.TrashFrame(self, datasets, annotationpath, gene=gene)
            viewWindow.Show()
        elif not datasets:
            transit_tools.show_error_dialog("Error: No datasets selected.")
            return
        else:
            transit_tools.show_error_dialog("Error: No annotation file selected.")
            return

    def scatterFunc(self, event):
        """ """
        # annotationpath = self.annotation
        datasets = self.ctrlSelected() + self.expSelected()
        if len(datasets) == 2:
            if self.verbose:
                transit_tools.log(
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
            transit_tools.show_error_dialog("Please make sure only two datasets are selected (across control and experimental datasets).")

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
            transit_tools.show_error_dialog(error_text)
            traceback.print_exc()
    

    def ExportSelectFunc(self, selected_name, test=""):
        # X = self.methodChoice.GetCurrentSelection()
        # selected_name = self.methodChoice.GetString(X)

        if self.verbose:
            transit_tools.log(
                "Selected Export Method: %s" % (selected_name)
            )

        for name in export_methods:
            if export_methods[name].label == selected_name:
                methodobj = export_methods[name].method
                try:
                    M = methodobj.from_gui(self)
                    if M:
                        thread = threading.Thread(target=M.Run())
                        thread.setDaemon(True)
                        thread.start()
                except Exception as e:
                    transit_tools.log("Error: %s" % str(e))
                    traceback.print_exc()

    def ConvertSelectFunc(self, selected_name, test=""):
        annotationpath = self.annotation

        for name in convert_methods:
            if convert_methods[name].label == selected_name:
                methodobj = convert_methods[name].method
                try:
                    M = methodobj.from_gui(self)
                    if M:
                        thread = threading.Thread(target=M.Run())
                        thread.setDaemon(True)
                        thread.start()
                except Exception as e:
                    transit_tools.log("Error: %s" % str(e))
                    traceback.print_exc()

    def LoessPrevFunc(self, event):
        datasets_selected = self.ctrlSelected() + self.expSelected()
        if not datasets_selected:
            transit_tools.show_error_dialog("Need to select at least one control or experimental dataset.")
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

    def annotationPT_to_GFF3(self, event):
        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".gff3"
        # defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()
        outputPath = self.SaveFile(defaultDir, defaultFile)

        ORGANISM = transit_tools.fetch_name(annotationpath)
        if not annotationpath:
            transit_tools.show_error_dialog("Error: No annotation file selected.")

        elif outputPath:
            if self.verbose:
                transit_tools.log(
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
                transit_tools.log("Finished conversion")

    def annotationPT_to_PTT(self, event):

        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".ptt.table"
        # defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        datasets = self.ctrlSelected() + self.expSelected()
        if not annotationpath:
            transit_tools.show_error_dialog("Error: No annotation file selected.")
        elif not datasets:
            transit_tools.show_error_dialog("Error: Please add a .wig dataset, to determine TA sites.")
        else:

            outputPath = self.SaveFile(defaultDir, defaultFile)
            if not outputPath:
                return
            if self.verbose:
                transit_tools.log(
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
                transit_tools.log("Finished conversion")

    def annotationPTT_to_PT(self, event):

        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".prot_table"
        # defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        datasets = self.ctrlSelected() + self.expSelected()
        if not annotationpath:
            transit_tools.show_error_dialog("Error: No annotation file selected.")
        # elif not datasets:
        #    transit_tools.show_error_dialog("Error: Please add a .wig dataset, to determine TA sites.")
        else:

            outputPath = self.SaveFile(defaultDir, defaultFile)
            if not outputPath:
                return
            if self.verbose:
                transit_tools.log(
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
                transit_tools.log("Finished conversion")

    def annotationGFF3_to_PT(self, event):

        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".prot_table"
        # defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        datasets = self.ctrlSelected() + self.expSelected()
        if not annotationpath:
            transit_tools.show_error_dialog("Error: No annotation file selected.")
        else:
            outputPath = self.SaveFile(defaultDir, defaultFile)
            if not outputPath:
                return
            if self.verbose:
                transit_tools.log(
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
                transit_tools.log("Finished conversion")