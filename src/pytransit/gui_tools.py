import traceback
from typing import NamedTuple

from pytransit.transit_tools import wx
from pytransit.core_data import universal



bit_map = None

def rgba(*args):
    return tuple(args)

def bind_to(wxPythonObj, event):
    """
        Usage:
            @bind_to(gui_object, wx.EVENT_THING)
            def callback(event):
                pass
        
    """
    def wrapper2(function_to_attach):
        wxPythonObj.Bind(event, function_to_attach)
        def wrapper1(*args, **kwargs):
            try:
               return function_to_attach(*args, **kwargs)
            except Exception as error:
                handle_error(error)
        return wrapper1
    return wrapper2

def handle_error(error_obj):
    """
        Summary:
            logs error message in bottom corner of GUI
            prints it to console
            and avoids crashing the whole runtime
        Usage:
            try:
                pass
            except Exception as error:
                handle_error
    """
    traceback.print_exc()
    frame = universal.frame
    if frame and hasattr(frame, "statusBar") and hasattr(frame.statusBar, "SetStatusText"):
        frame.statusBar.SetStatusText("Error: "+str(error_obj.args))

def handle_traceback(traceback_obj):
    import traceback
    frame = universal.frame
    print(''.join(traceback.format_tb(traceback_obj)))
    if frame and hasattr(frame, "statusBar") and hasattr(frame.statusBar, "SetStatusText"):
        frame.statusBar.SetStatusText("Error: "+str(error.args))

def show_message(MSG=""):
    # TODO: Write docstring
    wx.MessageBox(MSG, "Info", wx.OK | wx.ICON_INFORMATION)

def set_status(message):
    frame = universal.frame
    if frame:
        frame.statusBar.SetStatusText(message)

def ask_for_files(message):
    import os
    frame = universal.frame
    file_dialog = wx.FileDialog(
        frame,
        message=message,
        defaultDir=os.getcwd(),
        defaultFile="",
        style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
    )
    output = None
    if file_dialog.ShowModal() == wx.ID_OK:
        output = list(file_dialog.GetPaths())
    file_dialog.Destroy()
    return output

def ask_for_output_file_path(
        default_folder=None,
        default_file_name="",
        output_extensions=u'Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
    ):
        path = None
        if not default_folder:
            default_folder = os.getcwd()

        
        file_dialog = wx.FileDialog(
            self,
            message="Save file as ...",
            defaultDir=default_folder,
            defaultFile=default_file_name,
            wildcard=output_extensions,
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
        )
        if file_dialog.ShowModal() == wx.ID_OK:
            path = file_dialog.GetPath()
        file_dialog.Destroy()
        return path

class color(NamedTuple):
    black          =  rgba(0, 0, 0)
    white          =  rgba(255, 255, 255)
    light_gray     =  rgba(199, 203, 205)
    gray           =  rgba(84, 110, 122)
    rust           =  rgba(193, 126, 112)
    orange         =  rgba(247, 140, 108)
    yellow         =  rgba(254, 195, 85)
    bananna_yellow =  rgba(221, 215, 144)
    lime           =  rgba(195, 232, 141)
    green           = rgba(187, 237, 181, 255)
    bold_green     =  rgba(78, 201, 176, 209)
    vibrant_green  =  rgba(4, 216, 149)
    dim_green      =  rgba(128, 203, 171)
    dark_slate     =  rgba(63, 132, 141)
    light_slate    =  rgba(100, 186, 197)
    light_blue     =  rgba(137, 221, 255)
    blue           =  rgba(130, 170, 255)
    electric_blue  =  rgba(0, 174, 255, 232)
    light_purple   =  rgba(199, 146, 234)
    pink           =  rgba(229, 126, 179)
    red            =  rgba(255, 85, 114)
    soft_red       =  rgba(240, 113, 120)

class NiceErrorLog(object):
    """
        Example:
            with nice_error_log:
                print("stuff")
                raise Exception("example error")
                
                # ((error magically gets logged and caught))
    """
    def __init__(*args, **kwargs):
        pass
    
    def __enter__(self):
        frame = universal.frame
        return frame
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            print(''.join(traceback.format_tb(traceback_obj)))
            frame = universal.frame
            if frame:
                frame.statusBar.SetStatusText("Error: "+str(error.args))

nice_error_log = NiceErrorLog()