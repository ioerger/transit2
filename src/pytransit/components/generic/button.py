from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools

class Button:
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, text="", background_color=None, default_size=(250, -1), min_size=None, max_size=None):
        window       = gui_tools.window
        
        # 
        # wx_object
        # 
        wx_object = GenBitmapTextButton(
            window,
            1,
            gui_tools.bit_map,
            text,
            size=wx.Size(*default_size),
        )
        if background_color: wx_object.SetBackgroundColour(gui_tools.color.green)
        if max_size        : wx_object.SetMaxSize(wx.Size(*max_size))
        if min_size        : wx_object.SetMinSize(wx.Size(*min_size))
        
        self.wx_object = wx_object
        self.events = LazyDict(
            on_click=lambda func: wx_object.Bind(wx.EVT_BUTTON, func),
        )
        self._state = LazyDict(
            # None
        )
        
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)