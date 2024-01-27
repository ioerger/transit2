import os
import json

from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.generic_tools.misc import singleton, no_duplicates, flatten_once, human_readable_data, is_iterable
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename, working_directory
from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools

from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table

# 
# data shared between functions
# 
samples = LazyDict(
    sample_sizer=None,
    wig_dropdown_wxobj=None,
    wig_table=None,
    wig_dropdown_options={},
    condition_header_sizer=None,
    condition_dropdown_wxobj=None,
    conditions_table=None,
    condition_dropdown_options={},
    sample_button_creators = {
        # put "name": None here for buttons you want to be in a specific order
        # "LOESS": None,
    },
    load_button_args=dict(
        alignment=wx.ALIGN_LEFT,
        color=gui_tools.color.light_blue,
        size=(200, 40),
    ),
)

@gui.add_menu("File", "Clear Loaded Data")
def load_combined_wig_clear_function(event):
    samples.wig_table.clear()
    gui.combined_wigs.clear()
    samples.conditions_table.clear()
    
    samples.cwig_getter and samples.cwig_getter.set_label("")
    samples.metadata_getter and samples.metadata_getter.set_label("")
    samples.annotation_getter and samples.annotation_getter.set_label("")

# 
# 
# internal tools 
# 
# 

# 
# for updating the dropdown options
# 
def update_wig_area_dropdown(new_choices):
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

def update_condition_area_dropdown(new_choices):
    with gui_tools.nice_error_log:
        # hide the old one before showing the new one
        if samples.condition_dropdown_wxobj != None:
            samples.condition_dropdown_wxobj.Hide()
        
        # if not (None or empty)
        if new_choices: 
            new_choices = {
                "[Select Tool]": lambda event: None, # always have a none option, and always make it the first option
                **new_choices,
            }
            
            samples.condition_dropdown_wxobj = wx.Choice(
                gui.frame,
                wx.ID_ANY,
                wx.DefaultPosition,
                (120, -1),
                list(new_choices.keys()),
                0,
            )
            samples.condition_dropdown_wxobj.SetSelection(0)
            samples.condition_header_sizer.Add(samples.condition_dropdown_wxobj, proportion=0, flag=wx.EXPAND, border=gui_tools.default_padding)
            samples.condition_header_sizer.Layout()
            gui.frame.Layout()
            
            @gui_tools.bind_to(samples.condition_dropdown_wxobj, wx.EVT_CHOICE)
            def _(event):
                choice = samples.condition_dropdown_wxobj.GetString(samples.condition_dropdown_wxobj.GetCurrentSelection())
                # run the callback that corrisponds the the choice
                new_choices[choice](event)

# is only called in transit_gui.py
def create_sample_area(frame):
    from pytransit.components.panel_helpers import create_file_input
    
    wx_object = None
    with Column() as samples.sample_sizer:
        wx_object = samples.sample_sizer.wx_object
        samples.cwig_getter = None
        samples.metadata_getter = None
        samples.annotation_getter = None
        
        with Row() as button_container:
            with Column() as load_button_container:
                # 
                # cwig, metadata, and annotation button
                # 
                if True:
                    def update_values(*args):
                        load_combined_wigs_and_metadatas(
                            cwig_paths=[samples.cwig_getter()],
                            metadata_paths=[samples.metadata_getter()],
                            annotation_paths=[samples.annotation_getter()],
                        )
                            
                    samples.cwig_getter       = create_file_input(frame, load_button_container.wx_object, button_label="Select Combined Wig File", after_select=update_values, **samples.load_button_args, init_file_text="Must Select File",
                                                                  tooltip_text="Multiple .wig files are combined into a single combined_wig file using the **'transit export combined_wig'** command on the command line.")
                    samples.metadata_getter   = create_file_input(frame, load_button_container.wx_object, button_label="Select Metadata File",     after_select=update_values, **samples.load_button_args, init_file_text="Must Select File",
                                                                  tooltip_text="File describing each of the samples (as a spreadsheet, e.g. in Excel, which is then saved in tab-separated format).  Most commonly, there will be several replicates associated with each condition.")
                    samples.annotation_getter = create_file_input(frame, load_button_container.wx_object, button_label="Select Annotation File",   after_select=update_values, **samples.load_button_args, init_file_text="Must Select File",
                                                                   tooltip_text="File in 'prot_table' format: a tab-separated text file containing the information on each gene, such as coordinates, strand, ORF id, and gene name.")
            
            button_container.add(
                load_button_container,
                expand=True,
                proportion=4,
            )
        
                
            samples.sample_sizer.add(
                button_container,
                expand=True,
                proportion=0,
            )
        
        # padding
        samples.sample_sizer.wx_object.Add(10, 10)
        
        # 
        # text
        # 
        samples.sample_sizer.add(
            Text("Samples"),
            proportion=0,
        )
        
        samples.wig_header_sizer = wx.BoxSizer(wx.HORIZONTAL) 
        samples.sample_sizer.add(
            samples.wig_header_sizer,
            expand=True,
            proportion=0,
        )
        
        # 
        # samples.wig_table
        # 
        with Table(
            soft_size=True
        ) as samples.wig_table:
            # 
            # show wig-specific buttons
            # 
            @samples.wig_table.events.on_select
            def _(event):
                update_wig_area_dropdown(samples.wig_dropdown_options)
                for each in samples.sample_button_creators.values():
                    with gui_tools.nice_error_log:
                        each(samples.wig_table, samples.wig_header_sizer)
            
            samples.sample_sizer.add(
                samples.wig_table.wx_object,
                proportion=1, # 29 does something strange
                border=5,
                expand=True,
            )
        
        samples.sample_sizer.add(
            Text("Conditions"),
            proportion=0,
        )
        
        # 
        # condition_header
        #
        samples.condition_header_sizer = wx.BoxSizer(wx.HORIZONTAL) 
        samples.sample_sizer.add(
            samples.condition_header_sizer,
            expand=True,
            proportion=0,
        )
        
        # 
        # samples.conditions_table
        # 
        with Table(
            soft_size=True
        ) as samples.conditions_table:
        
            # 
            # show condition-specific buttons
            # 
            @samples.conditions_table.events.on_select
            def _(event):
                update_condition_area_dropdown(samples.condition_dropdown_options)
                for each in samples.sample_button_creators.values():
                    with gui_tools.nice_error_log:
                        each(samples.conditions_table, samples.condition_header_sizer)
                        
            
            samples.sample_sizer.add(
                samples.conditions_table.wx_object,
                proportion=1, # 29 does something strange
                border=5,
                expand=True,
                
            )
    
    samples.wig_header_sizer.Layout()
    gui.frame.Layout()
            
    # 
    # preload files if in debugging mode
    # 
    if debugging_enabled:
        from os import remove, getcwd
        load_combined_wigs_and_metadatas(
            [f"{root_folder}/src/pytransit/data/iron.transit/comwig.tsv"],
            [f"{root_folder}/src/pytransit/data/iron.transit/metadata.tsv"],
            [f"{root_folder}/src/pytransit/data/iron.transit/annotation.H37Rv.prot_table"],
        )
        
    return wx_object
    
def load_combined_wigs_and_metadatas(cwig_paths, metadata_paths, annotation_paths):
    """
        just a helper method
        once the "TODO: FOR DEBUGGING ONLY" is removed, this should 
        probably be inlined again for clarity (will only have one caller)
    """
    cwig_paths       = [ each for each in cwig_paths       if each != None ]
    metadata_paths   = [ each for each in metadata_paths   if each != None ]
    annotation_paths = [ each for each in annotation_paths if each != None ]
    with gui_tools.nice_error_log:
        if cwig_paths:
            samples.wig_table.clear()
            samples.conditions_table.clear()
            gui.combined_wigs.clear()
        
        # 
        # load the data from the files
        # 
        for each_cwig_path, each_metadata_path, each_annotation_path in zip(cwig_paths, metadata_paths, annotation_paths):
            logging.log(f"Loading '{os.path.basename(each_cwig_path)}' and '{os.path.basename(each_metadata_path)}'")
            with gui_tools.nice_error_log:
                gui.combined_wigs.append(
                    tnseq_tools.CombinedWig.load(
                        main_path=each_cwig_path,
                        metadata_path=each_metadata_path,
                        annotation_path=each_annotation_path,
                    )
                )
        
        logging.log(f"Done")

            
        # 
        # add graphical entries for each sample and condition
        # 
        if len(gui.combined_wigs):
            number_of_new_combined_wigs = len(cwig_paths)
            combined_wig = gui.combined_wigs[-1]
            if number_of_new_combined_wigs > 0:
                from pytransit.generic_tools import misc
                for each_sample in combined_wig.samples:
                    # BOOKMARK: here's where "density", "nz_mean", and "total count" can be added (they just need to be calculated)
                    samples.wig_table.add({
                        # add hidden link to object
                        "__wig_obj":each_sample,
                        # NOTE: all of these names are used by other parts of the code (caution when removing or renaming them)
                        "id":each_sample.id,
                        "conditions":(",".join(each_sample.condition_names) or "[None]"),
                        "density":round(each_sample.extra_data.density, 4),
                        "total_insertions":round(each_sample.extra_data.sum),
                        "non_zero_mean":round(each_sample.extra_data.non_zero_mean),
                        
                        # # uncomment to show covariate values next to each sample
                        # # add in multiple columns, one per covariate
                        # **combined_wig.metadata.covariates_dict_for(
                        #     wig_fingerprint=each_sample.fingerprint
                        # ),
                        
                        # # uncomment to add additional summary data
                        # non_zero_median=each_sample.extra_data.non_zero_median,
                        # count=each_sample.extra_data.count,
                        # mean=each_sample.extra_data.mean,
                        # max=each_sample.extra_data.max,
                        # skew=each_sample.extra_data.skew,
                        # kurtosis=each_sample.extra_data.kurtosis,
                    })
                
                covariate_headers = combined_wig.metadata.covariate_names
                for each_condition in combined_wig.conditions:
                    rows_with_condition = [ row for row in gui.combined_wigs[-1].metadata.rows if row["Condition"] == each_condition.name ]
                    
                    samples.conditions_table.add({
                        "Condition": each_condition.name,
                        **{
                            each_header: ", ".join(misc.no_duplicates([
                               f"{each_row[each_header]}" for each_row in rows_with_condition 
                            ]))
                                for each_header in covariate_headers
                        },
                    })

# should only be called/used by globals.py
def get_selected_samples():
    return [ each["__wig_obj"] for each in samples.wig_table.selected_rows ]

def get_selected_condition_names():
    selected_rows = samples.conditions_table.selected_rows
    return no_duplicates([ each["name"] for each in samples.conditions_table.selected_rows ])

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
        samples.wig_dropdown_options[name] = wrapper
        return wrapper
    return decorator

def add_condition_area_dropdown_option(name):
    """
    Example Usage:
        @gui.add_condition_area_dropdown_option(name="Click Me", size=(120, -1))
        def on_button_click(event):
            print("Howdy")
    """
    def decorator(function_being_wrapped):
        def wrapper(*args,**kwargs):
            with gui_tools.nice_error_log:
                return function_being_wrapped(*args, **kwargs)
        samples.condition_dropdown_options[name] = wrapper
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

