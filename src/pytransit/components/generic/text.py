from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename
from pytransit.specific_tools import  gui_tools

class Text:
    """
        Overview:
            self.wx_object
    """
    def __init__(self, content, font_size=12, underline=False, bold=False, italic=False):
        frame       = gui.frame
        wx_object = wx.StaticText( # not actually static btw
            frame,
            wx.ID_ANY,
            content,
            wx.DefaultPosition,
            wx.DefaultSize,
            0,
        )
        
        font_info = wx.FontInfo(font_size)
        if bold     : font_info = font_info.Bold()
        if underline: font_info = font_info.Underline()
        if italic   : font_info = font_info.Italic()
        # TODO: font family: https://docs.wxpython.org/wx.FontFamily.enumeration.html#wx-fontfamily
        wx_object.SetFont(wx.Font(font_info))
        
        self.wx_object = wx_object
        self.events = LazyDict(
            # None
        )
        self._state = LazyDict(
            content=content,
        )
    
    @property
    def content(self):
        return self._state.content
    
    @content.setter
    def content(self, value):
        self._state.content = value
        self.wx_object.SetLabel(value)
    
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)