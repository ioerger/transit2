import os

from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename, working_directory, read_result

import pytransit.components.file_display as file_display
from pytransit.specific_tools import  gui_tools

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table

results = LazyDict(
    header=None,
    table=None,
    file_action_choice_element=None,
    dropdown_width=220,
)
def create_results_area(frame):
    results_sizer = wx.BoxSizer(wx.VERTICAL)
    
    # 
    # Box
    # 
    if True:
        results_sizer.Add(10, 20) # spacer
        results_sizer.Add(
            wx.StaticText(
                gui.frame,
                label="Results Files",
                size=(330, -1),
            ),
            proportion=0,
            flag=wx.LEFT,
            border=gui_tools.default_padding,
        )
        
        results.header = wx.BoxSizer(wx.HORIZONTAL)
            
        # 
        # add_file_button
        # 
        if True:
            add_file_button = GenBitmapTextButton(
                frame,
                1,
                gui_tools.bit_map,
                "Add Results File",
                size=wx.Size(150, 30)
            )
            results.header.Add(add_file_button, 0, wx.ALL, gui_tools.default_padding)
            
            @gui_tools.bind_to(add_file_button, wx.EVT_BUTTON)
            def _(event):
                with gui_tools.nice_error_log:
                    dlg = wx.FileDialog(
                        frame,
                        message="Choose a file",
                        defaultDir=working_directory,
                        defaultFile="",
                        wildcard="Results Files (*.dat)|*.dat;|\nResults Files (*.txt)|*.txt;|\nAll files (*.*)|*.*",
                        style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                    )
                    paths = []
                    if dlg.ShowModal() == wx.ID_OK:
                        paths = dlg.GetPaths()
                    dlg.Destroy()
                    # ^close the menu before processing the files to be more responsive
                    
                    print("You chose the following Results file(s):")
                    for fullpath in paths:
                        add(fullpath)
            
        # 
        # fileActionChoice
        # 
        if True:
            # by default no options are available
            change_file_action_choices({})

            
        results_sizer.Add(results.header, 0, 0, gui_tools.default_padding)
    
    # 
    # results.table
    # 
    with Table(
        initial_columns=[ "name", "type", "path"],
        soft_size=True,
    ) as results.table:
    
        @results.table.events.on_select
        def _(event):
            with gui_tools.nice_error_log:
                row = results.table.rows[event.GetIndex()]
                dropdown_options_for_row = row.get("__dropdown_options", None)
                if isinstance(dropdown_options_for_row, dict):
                    # attach all their callbacks to the dropdown
                    change_file_action_choices(dropdown_options_for_row)
            
        results_sizer.Add(
            results.table.wx_object,
            1,
            wx.ALL | wx.EXPAND,
            gui_tools.default_padding,
        )
        
    return results_sizer


def change_file_action_choices(new_choices):
    with gui_tools.nice_error_log:
            
        # hide the old one before showing the new one
        if results.file_action_choice_element != None:
            results.file_action_choice_element.Hide()
        
        # if not (None or empty)
        if new_choices: 
            new_choices = {
                "[Select Tool]": lambda event: None, # always have a none option, and always make it the first option
                **new_choices,
            }
            
            results.file_action_choice_element = wx.Choice(
                gui.frame,
                wx.ID_ANY,
                wx.DefaultPosition,
                (results.dropdown_width, -1),
                list(new_choices.keys()),
                0,
            )
            results.file_action_choice_element.SetSelection(0)
            results.header.Add(results.file_action_choice_element, proportion=0, flag=wx.EXPAND, border=gui_tools.default_padding)
            results.header.Layout()
            gui.frame.Layout()
            
            @gui_tools.bind_to(results.file_action_choice_element, wx.EVT_CHOICE)
            def _(event):
                choice = results.file_action_choice_element.GetString(results.file_action_choice_element.GetCurrentSelection())
                # run the callback that corrisponds the the choice
                with gui_tools.nice_error_log:
                    new_choices[choice](event)
    

def add(path):
    if gui.is_active:
        with gui_tools.nice_error_log:
            # if not a recognized file type
            values_for_result_table = dict(
                name=basename(path),
                type="Unknown",
                path=path,
            )
            # if recognized
            result_object = read_result(path)
            if result_object:
                values_for_result_table = result_object.values_for_result_table
            
            if debugging_enabled:
                print(LazyDict(
                    result_object=str(result_object),
                    values_for_result_table=values_for_result_table,
                ))
            results.table.add(values_for_result_table)