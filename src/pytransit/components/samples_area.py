import os

from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename, working_directory
import pytransit.gui_tools as gui_tools

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table

from pytransit.components.comwig_picker import create_comwig_picker


# 
# Samples
# 
wig_table = None
def create_sample_area(frame):
    global wig_table
    
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
        # combined_wig_file_picker
        # 
        if True:
            # 
            # component
            # 
            combined_wig_file_picker = GenBitmapTextButton(
                frame,
                1,
                gui_tools.bit_map,
                "Add Files",
                size=wx.Size(250, -1),
            )
            combined_wig_file_picker.SetBackgroundColour(gui_tools.color.green)
            
            # 
            # callback
            # 
            @gui_tools.bind_to(combined_wig_file_picker, wx.EVT_BUTTON)
            def load_combined_wig_file_func(event): # BOOKMARK: cwig_callback
                file_dialog = wx.FileDialog(
                    frame,
                    message="Choose a cwig file",
                    defaultDir=working_directory,
                    defaultFile="",
                    wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                    style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                )
                if file_dialog.ShowModal() == wx.ID_OK:
                    cwig_paths = list(file_dialog.GetPaths())
                    metadata_paths = []
                    for fullpath in cwig_paths:
                        metadata_dialog = wx.FileDialog(
                            frame,
                            message=f"\n\nPick the sample metadata\nfor {basename(fullpath)}\n\n",
                            defaultDir=working_directory,
                            defaultFile="",
                            wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                            style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                        )
                        if metadata_dialog.ShowModal() == wx.ID_OK:
                            metadata_path = metadata_dialog.GetPaths()[0]
                            metadata_paths.append(
                                metadata_path
                            )
                        
                        metadata_dialog.Destroy()
                    for each_cwig_path, each_metadata_path in zip(cwig_paths, metadata_paths):
                        universal.session_data.add_cwig(
                            cwig_path=each_cwig_path,
                            metadata_path=each_metadata_path,
                        )
                file_dialog.Destroy()
                
                # 
                # add graphical entries for each condition
                # 
                if True:
                    for each_condition in universal.session_data.conditions:
                        wig_table.add(dict(
                            name=each_condition.name,
                            disabled=each_condition.is_disabled,
                        ))
            
            
            inner_sample_sizer.Add(combined_wig_file_picker, 1, wx.ALIGN_CENTER_VERTICAL, 5)
            
        outer_sample_sizer.Add(inner_sample_sizer, 0, wx.EXPAND, 5)
    
    # 
    # wig_table
    # 
    with Table() as wig_table:
        outer_sample_sizer.Add(
            wig_table.wx_object,
            1, # id?
            wx.ALL | wx.EXPAND,
            5 # ??
        )
    
    return outer_sample_sizer