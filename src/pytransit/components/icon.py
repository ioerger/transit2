from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton

from pytransit.globals import logging
from pytransit.specific_tools import gui_tools

if not HAS_WX:
    class InfoIcon: pass
else:
    class InfoIcon(wx.StaticBitmap):
        def __init__(self, panel, flag, bit_map=None, tooltip=""):
            if not bit_map:
                bit_map = wx.ArtProvider.GetBitmap(
                    wx.ART_INFORMATION, wx.ART_OTHER, (16, 16)
                )
            wx.StaticBitmap.__init__(self, panel, flag, bit_map)
            tp = wx.ToolTip(tooltip)
            self.SetToolTip(tp)