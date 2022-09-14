import sys

from pytransit.tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub

import traceback
import datetime
from pytransit.tools import logging, transit_tools
from pytransit.components.menu import convert_menu_item
from pytransit.components.icon import InfoIcon

class ConvertGUI:
    def __init__(self):
        self.wxobj = None
        self.menuitem = None
        self.LABELSIZE = (100, -1)
        self.WIDGETSIZE = (100, -1)

    def define_menu_item(self, wxobj, label):
        self.wxobj = wxobj
        self.menuitem = wx.MenuItem(
            convert_menu_item,
            wx.ID_ANY,
            label,
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )

class ConvertMethod:
    def __init__(
        self,
        short_name,
        long_name,
        description,
        label,
        annotation_path,
        output,
        wxobj=None,
    ):
        self.short_name = short_name
        self.long_name = long_name
        self.description = description
        self.label = label
        self.output = output
        self.annotation_path = annotation_path

        self.WX_VERSION = WX_VERSION
        self.wxobj = wxobj

    @classmethod
    def from_gui(self, wxobj):
        raise NotImplementedError

    @classmethod
    def from_args(self, args, kwargs):
        raise NotImplementedError

    def Run(self):
        raise NotImplementedError

    def finish(self):
        if self.wxobj:
            wx.CallAfter(pub.sendMessage, "finish", msg=self.short_name.lower())

class TransitConvert:
    def __init__(
        self, sn, ln, desc, lab, method_class=ConvertMethod, gui_class=ConvertGUI
    ):
        self.short_name  = sn
        self.long_name   = ln
        self.description = desc
        self.short_desc  = self.long_name
        self.long_desc   = desc
        self.full_name   = f"[{self.short_name}]  -  {self.short_desc}"
        self.label       = lab
        self.method      = method_class
        self.gui         = gui_class()
        self.transposons_text = ""

    def __str__(self):
        return f"""
            Convert Method:
                Short Name:  {self.short_name}
                Long Name:   {self.long_name}
                Description: {self.description}
                Method:      {self.method}
                GUI:         {self.gui}
        """.replace('\n            ','\n').strip()