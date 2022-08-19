import os

from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename, working_directory
import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table


# 
# Samples
# 
sample_table, conditions_table = None, None
def create_sample_area(frame):
    global sample_table
    global conditions_table
    
    wx_object = None
    with Column() as outer_sample_sizer:
        wx_object = outer_sample_sizer.wx_object
        
        outer_sample_sizer.add(
            Text("Samples"),
            proportion=0,
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
                    with gui_tools.nice_error_log:
                        if not universal.get("busy_running_method", False): # apparently this hook triggers for ALL button presses, so we must filter for when THIS button was clicked 
                            file_dialog = wx.FileDialog(
                                frame,
                                message="Choose a cwig file",
                                defaultDir=working_directory,
                                defaultFile="",
                                # wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                                style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                            )
                            cwig_paths = []
                            metadata_paths = []
                            if file_dialog.ShowModal() == wx.ID_OK:
                                cwig_paths = list(file_dialog.GetPaths())
                                metadata_paths = []
                                for fullpath in cwig_paths:
                                    metadata_dialog = wx.FileDialog(
                                        frame,
                                        message=f"\n\nPick the sample metadata\nfor {basename(fullpath)}\n\n",
                                        defaultDir=working_directory,
                                        defaultFile="",
                                        # wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                                        style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                                    )
                                    if metadata_dialog.ShowModal() == wx.ID_OK:
                                        metadata_path = metadata_dialog.GetPaths()[0]
                                        metadata_paths.append(
                                            metadata_path
                                        )
                                    
                                    metadata_dialog.Destroy()
                            file_dialog.Destroy()
                            
                            # 
                            # load the data from the files
                            # 
                            for each_cwig_path, each_metadata_path in zip(cwig_paths, metadata_paths):
                                transit_tools.log(f"Loading '{os.path.basename(each_cwig_path)}' and '{os.path.basename(each_metadata_path)}'")
                                with gui_tools.nice_error_log:
                                    universal.session_data.combined_wigs.append(
                                        tnseq_tools.CombinedWig(
                                            main_path=each_cwig_path,
                                            metadata_path=each_metadata_path,
                                        )
                                    )
                            
                            transit_tools.log(f"Done")

                                
                            # 
                            # add graphical entries for each condition
                            # 
                            if True:
                                for each_sample in universal.session_data.samples:
                                    sample_table.add(dict(
                                        # NOTE: all of these names are used by other parts of the code (caution when removing or renaming them)
                                        name=basename(each_sample.path),
                                        condition=each_sample.extra_data.get("condition", "[None]"),
                                        path=each_sample.path,
                                    ))
                                
                                for each_condition in universal.session_data.conditions:
                                    conditions_table.add(dict(
                                        name=each_condition.name,
                                    ))
                
                inner_sample_sizer.Add(combined_wig_file_picker, 1, wx.ALIGN_CENTER_VERTICAL, 5)
                
            outer_sample_sizer.add(
                inner_sample_sizer,
                expand=True,
                proportion=0,
            )
        
        # 
        # sample_table
        # 
        with Table() as sample_table:
            outer_sample_sizer.add(
                sample_table.wx_object,
                proportion=1, # 29 does something strange
                border=5,
                expand=True,
            )
        
        outer_sample_sizer.add(
            Text("Conditions"),
            proportion=0,
        )
        
        # 
        # conditions_table
        # 
        with Table() as conditions_table:
            outer_sample_sizer.add(
                conditions_table.wx_object,
                proportion=1, # 29 does something strange
                border=5,
                expand=True,
            )
        
    return wx_object