import os

from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools import gui_tools

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button

annotation_data = None
# 
# annotation_wrapper
# 
def create_annotation_area(frame):
    global annotation_data
    annotation_wrapper = wx.StaticBoxSizer(
        wx.StaticBox(
            frame,
            wx.ID_ANY,
            ""
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
                "Annotation File:",
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
                wildcard="prot_table or GFF3 files (*.gff3;*.gff;*.prot_table;*.txt)|*.gff3;*.gff;*.prot_table;*.txt",
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
            def annotation_file_func(event):
                global annotation_data
                annotation_data = event.GetPath()
                gui._annotation_path = annotation_data
                for each_combined_wig in gui.combined_wigs:
                    each_combined_wig.annotation_path = gui.annotation_path # TODO: this will probably change in the future to allow each combined_wig to have its own annotation_path

        annotation_wrapper.Add(annot_sizer, 1, wx.EXPAND, 5)
    
    return annotation_wrapper
    