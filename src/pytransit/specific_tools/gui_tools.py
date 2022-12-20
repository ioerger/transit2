import traceback
import os
from typing import NamedTuple

from pytransit.specific_tools.transit_tools import wx
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled



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
    frame = gui.frame
    error_message = ''.join(traceback.format_tb(traceback_obj))
    print(error_message)
    if frame and hasattr(frame, "status_bar") and hasattr(frame.status_bar, "SetStatusText"):
        frame.status_bar.SetStatusText(error_message)

def set_status(message):
    frame = gui.frame
    if frame and gui.is_active and hasattr(frame, "status_bar"):
        frame.status_bar.SetStatusText(message)
        wx.Yield()

def ask_for_files(
        message, 
        *,
        default_folder=None,
        default_file_name="",
        allowed_extensions='All files (*.*)|*.*',
    ):
        if isinstance(allowed_extensions, list):
            string = 'Common extensions '
            string += "(" + ",".join([ f"*.{each}" for each in allowed_extensions]) + ")|"
            string += "".join([ f"*.{each};" for each in allowed_extensions]) + "|"
            string += "\nAll files (*.*)|*.*"
            allowed_extensions = string
        
        import os
        output = []
        file_dialog = wx.FileDialog(
            gui.frame,
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
        if isinstance(allowed_extensions, list):
            string = 'Common extensions '
            string += "(" + ",".join([ f"*.{each}" for each in allowed_extensions]) + ")|"
            string += "".join([ f"*.{each};" for each in allowed_extensions]) + "|"
            string += "\nAll files (*.*)|*.*"
            allowed_extensions = string
        
        path = None
        file_dialog = wx.FileDialog(
            gui.frame,
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
        
        if isinstance(output_extensions, list):
            string = 'Common output extensions '
            string += "(" + ",".join([ f"*.{each}" for each in output_extensions]) + ")|"
            string += "".join([ f"*.{each};" for each in output_extensions]) + "|"
            string += "\nAll files (*.*)|*.*"
            output_extensions = string
        
        file_dialog = wx.FileDialog(
            gui.frame,
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
    light_gray     =  rgba(225, 225, 225)
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
        frame = gui.frame
        return frame
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            print(''.join(traceback.format_tb(traceback_obj)))
            error_message = " ".join([ f"{each}" for each in error.args])
            print(error_message)
            frame = gui.frame
            if frame and hasattr(frame, "status_bar"):
                frame.status_bar.SetStatusText("Error: "+error_message)

nice_error_log = NiceErrorLog()

def show_image(path):
    from pytransit.globals import logging, gui
    if gui.is_active:
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

# 
# Image Converters
# 
if True:
    def wx_bitmap_to_wx_image(my_bitmap):
        return wx.ImageFromBitmap(my_bitmap)

    def wx_bitmap_to_pil_image(my_bitmap):
        return wx_image_to_pil_image(wx_bitmap_to_wx_image(my_bitmap))

    def pil_image_to_wx_bitmap(my_pil_image):
        return wx_image_to_wx_bitmap(pil_image_to_wx_image(my_pil_image))

    def pil_image_to_wx_image(my_pil_image):
        my_wx_image = wx.EmptyImage(my_pil_image.size[0], my_pil_image.size[1])
        try:
            my_wx_image.SetData(my_pil_image.convert("RGB").tostring())
        except:
            my_wx_image.SetData(my_pil_image.convert("RGB").tobytes())
        return my_wx_image

    def wx_image_to_wx_bitmap(my_wx_image):
        return my_wx_image.ConvertToBitmap()
    
    def wx_image_to_pil_image(wx_image):
        raise Exception(f'''This function (wx_image_to_pil_image) was never implemented''')