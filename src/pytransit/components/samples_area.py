import os

from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.generic_tools.misc import singleton, no_duplicates, flatten_once, human_readable_data
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename, working_directory
from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools

from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table

# 
# Samples
# 
sample_table, conditions_table = None, None
sample_button_creators = []
def get_selected_samples():
    return [ each["__wig_obj"] for each in sample_table.selected_rows ]

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
                    "Load Combined Wig and Metadata",
                    size=wx.Size(-1, -1),
                    
                )
                combined_wig_file_picker.SetBackgroundColour(gui_tools.color.green)
                
                # 
                # callback
                # 
                @gui_tools.bind_to(combined_wig_file_picker, wx.EVT_BUTTON)
                def load_combined_wig_file_func(event):
                    with gui_tools.nice_error_log:
                        if not gui.busy_running_method: # apparently this hook triggers for ALL button presses, so we must filter for when THIS button was clicked 
                            file_dialog = wx.FileDialog(
                                frame,
                                message="Choose a cwig file",
                                defaultDir=working_directory,
                                defaultFile="",
                                # wildcard="Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
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
                                        # wildcard="Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                                        style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                                    )
                                    if metadata_dialog.ShowModal() == wx.ID_OK:
                                        metadata_path = metadata_dialog.GetPaths()[0]
                                        metadata_paths.append(
                                            metadata_path
                                        )
                                    
                                    metadata_dialog.Destroy()
                            file_dialog.Destroy()
                            
                            load_combined_wigs_and_metadatas(cwig_paths, metadata_paths)
                
                # inner_sample_sizer.Add(combined_wig_file_picker, proportion=2, flag=wx.ALIGN_LEFT, border=5)
                inner_sample_sizer.Add(combined_wig_file_picker)
                
            outer_sample_sizer.add(
                inner_sample_sizer,
                expand=True,
                proportion=0,
            )
        
        # 
        # sample_table
        # 
        with Table(
            max_size=(int(gui.width*0.7), 200)
        ) as sample_table:
            # 
            # show wig-specific buttons
            # 
            show_table_button = None
            @sample_table.events.on_select
            def _(event):
                for each in sample_button_creators:
                    with gui_tools.nice_error_log:
                        each(sample_table, inner_sample_sizer)
            
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
        with Table(
           max_size=(int(gui.width*0.7), 200)
        ) as conditions_table:
            outer_sample_sizer.add(
                conditions_table.wx_object,
                proportion=1, # 29 does something strange
                border=5,
                expand=True,
                
            )
    
    # 
    # preload files if in debugging mode
    # 
    if debugging_enabled:
        from os import remove, getcwd
        load_combined_wigs_and_metadatas(
            [f"{root_folder}/src/pytransit/data/111_cholesterol_glycerol_combined.cwig"],
            [f"{root_folder}/src/pytransit/data/222_samples_metadata_cg.txt"],
        )
        
    return wx_object
    
def load_combined_wigs_and_metadatas(cwig_paths, metadata_paths):
    """
        just a helper method
        once the "FIXME: FOR DEBUGGING ONLY" is removed, this should 
        probably be inlined again for clarity (will only have one caller)
    """
    
    # 
    # load the data from the files
    # 
    for each_cwig_path, each_metadata_path in zip(cwig_paths, metadata_paths):
        logging.log(f"Loading '{os.path.basename(each_cwig_path)}' and '{os.path.basename(each_metadata_path)}'")
        with gui_tools.nice_error_log:
            gui.combined_wigs.append(
                tnseq_tools.CombinedWig(
                    main_path=each_cwig_path,
                    metadata_path=each_metadata_path,
                )
            )
    
    logging.log(f"Done")

        
    # 
    # add graphical entries for each condition
    # 
    if True:
        for each_sample in gui.samples:
            # BOOKMARK: here's where "density", "nz_mean", and "total count" can be added (they just need to be calculated)
            sample_table.add(dict(
                # add hidden link to object
                __wig_obj=each_sample,
                # NOTE: all of these names are used by other parts of the code (caution when removing or renaming them)
                id=each_sample.id,
                conditions=(",".join(each_sample.condition_names) or "[None]"),
                density=round(each_sample.extra_data.density, 4),
                total_insertions=round(each_sample.extra_data.sum),
                non_zero_mean=round(each_sample.extra_data.non_zero_mean),
                # # uncomment to add additional summary data
                # non_zero_median=each_sample.extra_data.non_zero_median,
                # count=each_sample.extra_data.count,
                # mean=each_sample.extra_data.mean,
                # max=each_sample.extra_data.max,
                # skew=each_sample.extra_data.skew,
                # kurtosis=each_sample.extra_data.kurtosis,
            ))
        
        for each_condition in gui.conditions:
            conditions_table.add(dict(
                name=each_condition.name,
            ))


# 
# Let methods add Buttons
# 
def create_sample_area_button(name, size=(150, -1)):
    show_thing_button = None
    def decorator(function_being_wrapped):
        def create_show_thing_button(sample_table, inner_sample_sizer):
            nonlocal show_thing_button
            with gui_tools.nice_error_log:
                # hide an old button if it exists
                if show_thing_button != None:
                    show_thing_button.Hide()
                
                show_thing_button = GenBitmapTextButton(
                    gui.frame,
                    1,
                    gui_tools.bit_map,
                    name,
                    size=wx.Size(*size),
                )
                show_thing_button.SetBackgroundColour(gui_tools.color.light_blue)
                
                # 
                # callback
                # 
                @gui_tools.bind_to(show_thing_button, wx.EVT_BUTTON)
                def click_show_thing(event):
                    return function_being_wrapped(event)
            
            inner_sample_sizer.Add(show_thing_button)
            inner_sample_sizer.Layout()
        
        sample_button_creators.append(create_show_thing_button)
        return create_show_thing_button
    return decorator
