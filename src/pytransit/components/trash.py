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

try: import wx
except: pass

import pytransit.components.view_trash as view_trash
import pytransit.components.draw_trash as draw_trash
from pytransit.globals import logging
from math import *
import os
import ntpath
import traceback
from PIL import Image, ImageDraw, ImageFont

from pytransit.specific_tools import transit_tools, tnseq_tools, norm_tools, gui_tools

track_prefix = "[TrackView]"

# inherit from the MainFrame created in wxFowmBuilder and create CalcFrame
class TrashFrame(view_trash.MainFrame):
    # constructor
    def __init__(
        self,
        parent,
        wig_ids=["H37Rv_Sassetti_glycerol.wig"], # TODO: probably should remove these (I didnt add them) --Jeff
        annotation="H37Rv.prot_table",           # TODO: probably should remove these (I didnt add them) --Jeff
        gene="",
        scale=None,
        feature_hashes=[],
        feature_data=[],
    ):

        view_trash.MainFrame.__init__(self, parent)

        self.parent = parent
        self.size = wx.Size(1500, 800)
        self.start = 1
        self.end = 10000

        self.orf2data = tnseq_tools.AnnotationFile(path=annotation).orf_to_info
        self.hash = transit_tools.get_pos_hash(annotation)

        self.features = []

        # Data to facilitate search
        self.name2id = {}
        for orf, (name, desc, start, end, strand) in self.orf2data.items():
            name = name.lower()
            if name not in self.name2id:
                self.name2id[name] = []
            self.name2id[name].append(orf)

        self.lowerid2id = dict([(x.lower(), x) for x in self.orf2data.keys()])
        self.labels = wig_ids + ["All"]
        
        from pytransit.specific_tools.transit_tools import gather_sample_data_for
        self.fulldata, self.position = gather_sample_data_for(selected_samples=True)

        # Save normalized data
        (self.fulldata_norm, self.factors) = norm_tools.normalize_data(
            self.fulldata, method="nzmean"
        )
        self.wasNorm = False

        # initialize parent class

        self.feature_hashes = feature_hashes
        self.feature_data = feature_data

        if not scale:
            scale = [150] * len(wig_ids)
        self.scale = scale
        self.globalScale = False

        self.datasetChoice.SetItems(self.labels)
        self.datasetChoice.SetSelection(0)

        if gene:
            self.searchText.SetValue(gene)
            self.searchFunc(gene)

        self.updateFunc(parent)
        self.Fit()
        self.datasetChoice.SetSelection(len(self.labels) - 1)

    def updateFunc(self, event):
        try:
            self.DrawCanvas()
        except Exception as e:
            logging.log("ERROR: %s" % e)
            traceback.print_exc()

    def leftFunc(self, event):
        start = int(self.startText.GetValue())
        end = int(self.endText.GetValue())
        delta = end - start
        new_start = int(start - delta * 0.10)
        new_end = int(end - delta * 0.10)

        self.startText.SetValue(str(new_start))
        self.endText.SetValue(str(new_end))
        self.end = new_end
        self.start = new_start
        self.updateFunc(event)

    def rightFunc(self, event):
        start = int(self.startText.GetValue())
        end = int(self.endText.GetValue())
        delta = end - start
        new_start = int(start + delta * 0.10)
        new_end = int(end + delta * 0.10)
        self.end = new_end
        self.start = new_start

        self.startText.SetValue(str(new_start))
        self.endText.SetValue(str(new_end))
        self.updateFunc(event)

    def zoomInFunc(self, event):
        start = int(self.startText.GetValue())
        end = int(self.endText.GetValue())
        delta = end - start
        new_start = int(start + delta * 0.10)
        new_end = int(end - delta * 0.10)

        self.startText.SetValue(str(new_start))
        self.endText.SetValue(str(new_end))
        self.end = new_end
        self.start = new_start
        self.updateFunc(event)

    def zoomOutFunc(self, event):
        start = int(self.startText.GetValue())
        end = int(self.endText.GetValue())
        delta = end - start
        new_start = int(start - delta * 0.10)
        new_end = int(end + delta * 0.10)

        self.startText.SetValue(str(new_start))
        self.endText.SetValue(str(new_end))
        self.end = new_end
        self.start = new_start
        self.updateFunc(event)

    def changedMaxFunc(self, event):
        start = int(self.startText.GetValue())
        if self.maxText.GetValue() and self.maxText.GetValue() != "-":
            dataset_ii = self.datasetChoice.GetCurrentSelection()
            if dataset_ii == (len(self.labels) - 1):
                maxVal = int(self.maxText.GetValue())
                self.scale = [maxVal for _ in self.scale]
                logging.log("All Datasets scaled to %s" % (maxVal))
            else:
                self.scale[dataset_ii] = int(self.maxText.GetValue())
                logging.log(
                    "Dataset '%s' scaled to %s"
                    % (self.datasetChoice.GetString(dataset_ii), self.scale[dataset_ii])
                )
        self.updateFunc(event)

    def scaleFunc(self, event):
        if self.autoScaleCheck.GetValue():
            self.maxText.Enable(False)
            self.datasetChoice.Enable(False)
            self.globalScale = True
            logging.log("Scaling read-counts to local (Window) Maximum.")
        else:
            self.maxText.Enable(True)
            self.datasetChoice.Enable(True)
            self.globalScale = False
            logging.log("Scaling read-counts tracks individually.")
        self.updateFunc(event)

    def datasetSelectFunc(self, event):
        dataset_ii = self.datasetChoice.GetCurrentSelection()
        self.maxText.SetValue(str(self.scale[dataset_ii]))

    def resetFunc(self, event):
        self.startText.SetValue("1")
        self.endText.SetValue("10000")
        self.scale = [150] * len(self.labels)
        self.datasetChoice.SetSelection(len(self.labels) - 1)
        self.maxText.SetValue(str(self.scale[0]))

        self.updateFunc(event)

    def saveImageFunc(self, event):
        finished_image = self.m_bitmap1.GetBitmap()
        output_path = self.SaveFile(DIR=".", FILE="reads_canvas.png")
        if output_path:
            finished_image.SaveFile(output_path, wx.BITMAP_TYPE_PNG)
            logging.log("Image saved to the following path: %s" % output_path)

    def addFeatureFunc(self, event):
        wc = "Known Annotation Formats (*.prot_table,*.gff3,*.gff)|*.prot_table;*.gff3;*.gff;|\nProt Table (*.prot_table)|*.prot_table;|\nGFF3 (*.gff,*.gff3)|*.gff;*.gff3;|\nAll files (*.*)|*.*"
        path = self.OpenFile(DIR=self.parent.workdir, FILE="", WC=wc)
        try:
            if path:
                if self.checkHMMFeature(path):
                    H, S = self.getHMMHash(path)
                    self.feature_hashes.append(H)
                    self.feature_data.append(S)
                else:
                    self.feature_hashes.append(transit_tools.get_pos_hash(path))
                    self.feature_data.append(tnseq_tools.AnnotationFile(path=path).orf_to_info)
                self.updateFunc(self.parent)
                self.Fit()
            else:
                logging.log("No feature added")
        except Exception as e:
            logging.log("ERROR: %s" % e)
            traceback.print_exc()

    def checkHMMFeature(self, path):
        try:
            if open(path).readline().startswith("#HMM - Sites"):
                return True
        except Exception as e:
            logging.log("ERROR: %s" % e)
            return False
        return False

    def getHMMHash(self, path):
        hash = {}
        fi = open(path)
        line = fi.readline()
        tmp = line.split("\t")
        # print(tmp)
        last_state = ""
        start = 0
        end = 0
        counts = {"ES": 0, "NE": 0, "GD": 0, "GA": 0}
        states2info = {}
        while line:
            if line.startswith("#"):
                line = fi.readline()
                tmp = line.split("\t")
                continue

            current_pos, current_state = tmp[0], tmp[-2]
            if last_state != current_state:
                if last_state:
                    state_name = "%s-%d" % (last_state, counts[last_state])
                    counts[last_state] += 1
                    states2info[state_name] = (last_state, last_state, start, end, "+")
                    for pos in range(start, end + 1):
                        if pos not in hash:
                            hash[pos] = []
                        hash[pos].append(state_name)
                        start = int(current_pos)
                last_state = current_state
            else:
                end = int(current_pos)

            line = fi.readline()
            tmp = line.split("\t")
        return hash, states2info

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
            logging.log("You added the following file: %s" % path)
        dlg.Destroy()
        return path

    def searchFunc(self, event):

        query = self.searchText.GetValue()

        genes_match_name = self.name2id.get(query.lower(), [])
        gene_match_orf = self.lowerid2id.get(query.lower(), None)
        gene_match_orf_w_c = self.lowerid2id.get(query.lower() + "c", None)

        combined_match = genes_match_name
        if gene_match_orf:
            combined_match += [gene_match_orf]
        if gene_match_orf_w_c:
            combined_match += [gene_match_orf_w_c]

        if combined_match:
            logging.log("Genes matching query: %s" % ", ".join(combined_match))
        else:
            logging.log("No genes matching query!")

        if len(genes_match_name) == 1:  # Check if query is a name
            orf_match = genes_match_name[0]
        elif gene_match_orf:  # Check if query is an orf id
            orf_match = gene_match_orf
        elif gene_match_orf_w_c:  # Check if query is an orf id with missing "c" at end
            orf_match = gene_match_orf_w_c
        else:
            return

        (name, desc, start, end, strand) = self.orf2data.get(
            orf_match, ["-", "-", 0, 2000, "+"]
        )

        try:

            # Get min/max read info from text controls
            self.startText.SetValue(str(start))
            self.endText.SetValue(str(end))

            # Set the start and end coords with the new info
            self.start = start
            self.end = end
            self.DrawCanvas()

        except Exception as e:
            logging.log("ERROR: %s" % e)
            traceback.print_exc()

    def DrawCanvas(self):

        # self.autoScale = self.autoScaleCheck.GetValue()
        if self.normCheck.GetValue():
            image_pil = draw_trash.draw_canvas(
                self.fulldata_norm,
                self.position,
                self.hash,
                self.orf2data,
                self.feature_hashes,
                self.feature_data,
                labels=self.labels,
                scale=self.scale,
                globalScale=self.globalScale,
                start=self.start,
                end=self.end,
            )
            if not self.wasNorm:
                logging.log("Normalization factors: %s" % self.factors.flatten())
            self.wasNorm = True
        else:
            image_pil = draw_trash.draw_canvas(
                self.fulldata,
                self.position,
                self.hash,
                self.orf2data,
                self.feature_hashes,
                self.feature_data,
                labels=self.labels,
                scale=self.scale,
                globalScale=self.globalScale,
                start=self.start,
                end=self.end,
            )
            self.wasNorm = False

        image_wxImg = gui_tools.pil_image_to_wx_image(image_pil)
        self.m_bitmap1.SetBitmap(wx.BitmapFromImage(image_wxImg))
        self.Refresh()
        image_pil = ""

    def SaveFile(self, DIR=".", FILE="", WC=""):
        """
        Create and show the Save FileDialog
        """
        path = ""
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
        dlg.Destroy()
        return path


if __name__ == "__main__":
    app = wx.App(False)

    frame = TrashFrame(None)
    frame.Show(True)
    app.MainLoop()
