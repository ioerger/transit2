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
# data shared between functions
# 
samples = LazyDict(
    wig_header_sizer=None,
    wig_dropdown_wxobj=None,
    wig_table=None,
    dropdown_options={},
    conditions_table=None,
    sample_button_creators = {
        # put "name": None here for buttons you want to be in a specific order
        # "LOESS": None,
    },
)

# 
# 
# internal tools 
# 
# 

# 
# for updating the dropdown options
# 
def update_sample_area_dropdown(new_choices):
    with gui_tools.nice_error_log:
        # hide the old one before showing the new one
        if samples.wig_dropdown_wxobj != None:
            samples.wig_dropdown_wxobj.Hide()
        
        # if not (None or empty)
        if new_choices: 
            new_choices = {
                "[Select Tool]": lambda event: None, # always have a none option, and always make it the first option
                **new_choices,
            }
            
            samples.wig_dropdown_wxobj = wx.Choice(
                gui.frame,
                wx.ID_ANY,
                wx.DefaultPosition,
                (120, -1),
                list(new_choices.keys()),
                0,
            )
            samples.wig_dropdown_wxobj.SetSelection(0)
            samples.wig_header_sizer.Add(samples.wig_dropdown_wxobj, proportion=0, flag=wx.EXPAND, border=gui_tools.default_padding)
            samples.wig_header_sizer.Layout()
            gui.frame.Layout()
            
            @gui_tools.bind_to(samples.wig_dropdown_wxobj, wx.EVT_CHOICE)
            def _(event):
                choice = samples.wig_dropdown_wxobj.GetString(samples.wig_dropdown_wxobj.GetCurrentSelection())
                # run the callback that corrisponds the the choice
                new_choices[choice](event)

# is only called in transit_gui.py
def create_sample_area(frame):
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
            samples.wig_header_sizer = wx.BoxSizer(wx.HORIZONTAL)
            
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
                    size=wx.Size(300, 100),
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
                            
                            print(f"Loading paths:{cwig_paths}, {metadata_paths}")
                            load_combined_wigs_and_metadatas(cwig_paths, metadata_paths)
                
                samples.wig_header_sizer.Add(combined_wig_file_picker, proportion=1, flag=wx.ALIGN_CENTER_VERTICAL, border=5)
                # samples.wig_header_sizer.Add(combined_wig_file_picker)
                # inner_sample_sizer.Add(combined_wig_file_picker, 1, wx.ALIGN_CENTER_VERTICAL, 5)
                
            outer_sample_sizer.add(
                samples.wig_header_sizer,
                expand=True,
                proportion=0,
            )
        
        # 
        # samples.wig_table
        # 
        with Table(
            max_size=(int(gui.width*0.7), 200)
        ) as samples.wig_table:
            # 
            # show wig-specific buttons
            # 
            show_table_button = None
            @samples.wig_table.events.on_select
            def _(event):
                update_sample_area_dropdown(samples.dropdown_options)
                for each in samples.sample_button_creators.values():
                    with gui_tools.nice_error_log:
                        each(samples.wig_table, samples.wig_header_sizer)
            
            outer_sample_sizer.add(
                samples.wig_table.wx_object,
                proportion=1, # 29 does something strange
                border=5,
                expand=True,
            )
        
        outer_sample_sizer.add(
            Text("Conditions"),
            proportion=0,
            
        )
        
        # 
        # samples.conditions_table
        # 
        with Table(
           max_size=(int(gui.width*0.7), 200)
        ) as samples.conditions_table:
            outer_sample_sizer.add(
                samples.conditions_table.wx_object,
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
        number_of_new_combined_wigs = len(cwig_paths)
        if number_of_new_combined_wigs > 0:
            from pytransit.generic_tools import misc
            new_samples = misc.flatten_once([ each.samples for each in gui.combined_wigs[-number_of_new_combined_wigs:] ])
            for each_sample in new_samples:
                # BOOKMARK: here's where "density", "nz_mean", and "total count" can be added (they just need to be calculated)
                samples.wig_table.add(dict(
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
            
            new_conditions = misc.flatten_once([ each.conditions for each in gui.combined_wigs[-number_of_new_combined_wigs:] ])
            for each_condition in new_conditions:
                samples.conditions_table.add(dict(
                    name=each_condition.name,
                ))

# should only be called/used by globals.py
def get_selected_samples():
    return [ each["__wig_obj"] for each in samples.wig_table.selected_rows ]

def get_selected_condition_names():
    return [ each["name"] for each in samples.conditions_table.selected_rows ]

# 
# 
# External tools
# 
# 

# 
# Let methods add dropdown options
# 
def add_wig_area_dropdown_option(name):
    """
    Example Usage:
        @gui.add_wig_area_dropdown_option(name="Click Me", size=(120, -1))
        def on_button_click(event):
            print("Howdy")
    """
    def decorator(function_being_wrapped):
        def wrapper(*args,**kwargs):
            with gui_tools.nice_error_log:
                return function_being_wrapped(*args, **kwargs)
        samples.dropdown_options[name] = wrapper
        return wrapper
    return decorator

# 
# Let methods add Buttons
# 
def add_wig_area_button(name, size=(150, -1)):
    """
    Example Usage:
        @samples_area.add_wig_area_button(name="Click Me", size=(120, -1))
        def on_button_click(event):
            print("Howdy")
    """
    show_thing_button = None
    def decorator(function_being_wrapped):
        def create_show_thing_button(wig_table, wig_header_sizer):
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
                    with gui_tools.nice_error_log:
                        return function_being_wrapped(event)
            
            samples.wig_header_sizer.Add(show_thing_button)
            samples.wig_header_sizer.Layout()
        
        samples.sample_button_creators[name] = create_show_thing_button
        return create_show_thing_button
    return decorator

