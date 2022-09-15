import traceback
import os
from typing import NamedTuple

from pytransit.tools.transit_tools import wx
from pytransit.universal_data import universal



bit_map = None
default_padding = 5


def rgba(*args):
    return tuple(args)

def bind_to(wx_python_obj, event, *args, **kwargs):
    """
        Usage:
            @bind_to(gui_object, wx.EVENT_THING)
            def callback(event):
                pass
        
    """
    def wrapper2(function_to_attach):
        wx_python_obj.Bind(event, function_to_attach, *args, **kwargs)
        def wrapper1(*inner_args, **inner_kwargs):
            with nice_error_log:
               return function_to_attach(*inner_args, **inner_kwargs)
        return wrapper1
    return wrapper2

def handle_traceback(traceback_obj):
    import traceback
    frame = universal.frame
    print(''.join(traceback.format_tb(traceback_obj)))
    if frame and hasattr(frame, "status_bar") and hasattr(frame.status_bar, "SetStatusText"):
        frame.status_bar.SetStatusText("Error: "+str(error.args))

def set_status(message):
    frame = universal.frame
    if frame and universal.interface == "gui" and hasattr(frame, "status_bar"):
        frame.status_bar.SetStatusText(message)
        wx.Yield()

def ask_for_files(
        message, 
        *,
        default_folder=None,
        default_file_name="",
        allowed_extensions='All files (*.*)|*.*',
    ):
        import os
        output = []
        file_dialog = wx.FileDialog(
            universal.frame,
            message=message,
            defaultDir=default_folder or os.getcwd(),
            defaultFile=default_file_name,
            wildcard=allowed_extensions,
            style=wx.FD_MULTIPLE | wx.FD_SHOW_HIDDEN | wx.FD_OPEN,
        )
        if file_dialog.ShowModal() == wx.ID_OK:
            output = list(file_dialog.GetPaths())
        file_dialog.Destroy()
        return output

def ask_for_file(
        message, 
        *,
        default_folder=None,
        default_file_name="",
        allowed_extensions='All files (*.*)|*.*',
    ):
        
        path = None
        file_dialog = wx.FileDialog(
            universal.frame,
            message=message,
            defaultDir=default_folder or os.getcwd(),
            defaultFile=default_file_name,
            wildcard=allowed_extensions,
            style=wx.FD_SHOW_HIDDEN | wx.FD_OPEN,
        )
        if file_dialog.ShowModal() == wx.ID_OK:
            path = file_dialog.GetPath()
        file_dialog.Destroy()
        return path

def ask_for_output_file_path(
        default_folder=None,
        default_file_name="",
        output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
    ):
        path = None
        if not default_folder:
            default_folder = os.getcwd()

        
        file_dialog = wx.FileDialog(
            universal.frame,
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
            if frame and hasattr(frame, "status_bar"):
                frame.status_bar.SetStatusText("Error: "+str(error.args))

nice_error_log = NiceErrorLog()

def show_image(path):
    class ImgFrame(wx.Frame):
        def __init__(self, parent, filePath):
            wx.Frame.__init__(
                self,
                parent,
                id=wx.ID_ANY,
                title="%s" % (filePath),
                pos=wx.DefaultPosition,
                size=wx.Size(1150, 740),
                style=wx.DEFAULT_FRAME_STYLE | wx.TAB_TRAVERSAL,
            )
            self.SetSizeHintsSz(wx.DefaultSize, wx.DefaultSize)
            bSizer1 = wx.BoxSizer(wx.VERTICAL)

            self.m_bitmap1 = wx.StaticBitmap(
                self, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0
            )
            bSizer1.Add(self.m_bitmap1, 1, wx.ALL | wx.EXPAND, 5)
            self.SetSizer(bSizer1)
            self.Layout()
            self.Centre(wx.BOTH)

            img = wx.Image(filePath, wx.BITMAP_TYPE_ANY)
            self.m_bitmap1.SetBitmap(wx.BitmapFromImage(img))

            self.Refresh()
            self.Fit()
    
    ImgFrame(None, path).Show()


def run_method_by_label(*, method_options, method_label):
    with transit_tools.nice_error_log:
        for name in method_options:
            method_option = method_options[name]
            if method_option.label == method_label:
                method_object = method_option.method.from_gui(frame)
                if method_object:
                    thread = threading.Thread(target=method_object.Run())
                    thread.setDaemon(True)
                    thread.start()
