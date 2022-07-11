from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools
from pytransit.components.generic.box import Column

class WindowManager:
    """
        Overview:
            self.wx_object
            self.add(python_obj)
    """
    def __init__(self, default_size=(-1, -1), min_size=(700, -1), scroll_rate=(5,5)):
        wx_object = wx.ScrolledWindow(
            self,
            wx.ID_ANY,
            wx.DefaultPosition,
            wx.Size(*size),
            wx.HSCROLL | wx.VSCROLL,
        )
        wx_object.SetScrollRate(*scroll_rate)
        wx_object.SetMinSize(wx.Size(*min_size))
        
        self.wx_object = wx_object
        self.events = LazyDict(
            # None
        )
        self._state = LazyDict(
            sizer = Column(),
        )
    
    def add(self, *args, **kwargs):
        self.sizer.add(*args, **kwargs)
        return self
    
    @property
    def length(self):
        return self.sizer.length
    
    def __len__(self):
        return self.length
    
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)