from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools

class Text:
    """
        Overview:
            self.wx_object
    """
    def __init__(self, content):
        window       = gui_tools.window
        self.wx_object = wx.StaticText( # not actually static btw
            self.mainWindow,
            wx.ID_ANY,
            content,
            wx.DefaultPosition,
            wx.DefaultSize,
            0,
        )
        
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