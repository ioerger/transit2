import os

from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table

from pytransit.components.comwig_picker import create_comwig_picker


# 
# Samples
# 
def create_sample_area(frame):
    outer_sample_sizer = wx.StaticBoxSizer(
        wx.StaticBox(
            frame,
            wx.ID_ANY,
            u"Samples"
        ),
        wx.VERTICAL,
    )
    
    # 
    # box
    # 
    if True:
        inner_sample_sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # 
        # file_picker
        # 
        if True:
            file_picker = create_comwig_picker(frame)
            inner_sample_sizer.Add(file_picker, 1, wx.ALIGN_CENTER_VERTICAL, 5)
            
        outer_sample_sizer.Add(inner_sample_sizer, 0, wx.EXPAND, 5)
    
    # 
    # wig_table
    # 
    with Table() as component:
        outer_sample_sizer.Add(
            component.wx_object,
            1, # id?
            wx.ALL | wx.EXPAND,
            5 # ??
        )
    
    return outer_sample_sizer