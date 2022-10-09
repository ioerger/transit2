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

import ntpath
import subprocess
import os
import sys
from functools import partial

from pytransit.specific_tools.transit_tools import wx

if wx:
    class ImgFrame(wx.Frame):
        def __init__(self, parent, filePath):
            wx.Frame.__init__(
                self,
                parent,
                id=wx.ID_ANY,
                title="%s" % (filePath),
                pos=wx.DefaultPosition,
                size=wx.Size(1150, 740),
                style=wx.DEFAULT_FRAME_STYLE | wx.TAB_TRAVERSAL,
            )
            self.SetSizeHintsSz(wx.DefaultSize, wx.DefaultSize)
            bSizer1 = wx.BoxSizer(wx.VERTICAL)

            self.m_bitmap1 = wx.StaticBitmap(
                self, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0
            )
            bSizer1.Add(self.m_bitmap1, 1, wx.ALL | wx.EXPAND, 5)
            self.SetSizer(bSizer1)
            self.Layout()
            self.Centre(wx.BOTH)

            img = wx.Image(filePath, wx.BITMAP_TYPE_ANY)
            self.m_bitmap1.SetBitmap(wx.BitmapFromImage(img))

            self.Refresh()
            self.Fit()
