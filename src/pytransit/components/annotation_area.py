import os

from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
from pytransit.core_data import SessionData, universal
import pytransit.gui_tools as gui_tools

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button




# 
# annotation_wrapper
# 
def create_annotation_area(frame):
    annotation_wrapper = wx.StaticBoxSizer(
        wx.StaticBox(
            frame,
            wx.ID_ANY,
            u""
        ),
        wx.VERTICAL,
    )
    
    # 
    # annotation sizer
    # 
    if True:
        annot_sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # 
        # text
        # 
        if True:
            label_annot = wx.StaticText(
                frame,
                wx.ID_ANY,
                u"Annotation File:",
                wx.DefaultPosition,
                wx.DefaultSize,
                0,
            )
            annot_sizer.Add(
                label_annot,
                0,
                wx.ALIGN_CENTER_VERTICAL,
                0
            )
        
        # 
        # picker
        # 
        if True:
            annotation_file_picker = wx.FilePickerCtrl(
                frame,
                id=wx.ID_ANY,
                size=(400, 30),
                wildcard=u"prot_table or GFF3 files (*.gff3;*.gff;*.prot_table;*.txt)|*.gff3;*.gff;*.prot_table;*.txt",
                message="Select Annotation file (.prot_table or .gff3)",
                style=wx.FLP_DEFAULT_STYLE | wx.FLP_USE_TEXTCTRL | wx.FD_MULTIPLE,
            )
            annotation_file_picker.SetInitialDirectory(os.getcwd())
            annot_sizer.Add(
                annotation_file_picker,
                proportion=1,
                flag=wx.EXPAND | wx.ALL,
                border=5,
            )
            
            @gui_tools.bind_to(annotation_file_picker, wx.EVT_FILEPICKER_CHANGED)
            def annotationFileFunc(event):
                frame.annotation = event.GetPath() # this is for reducing breakages
                universal.annotation = frame.annotation

        annotation_wrapper.Add(annot_sizer, 1, wx.EXPAND, 5)
    
    return annotation_wrapper
    