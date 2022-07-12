from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools
from pytransit.components.generic.box import Column, Row

class Frame:
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, parent, title, default_size=(1350, 975), min_size=None, max_size=None):
        window = gui_tools.window
        
        # 
        # wx_object
        # 
        wx.Frame.__init__(
            window,
            parent,
            id=wx.ID_ANY,
            title=title,
            pos=wx.DefaultPosition,
            size=wx.Size(*default_size),
            style=wx.DEFAULT_FRAME_STYLE | wx.TAB_TRAVERSAL,
        )
        
        self.events = LazyDict(
            # None
        )
        self._state = LazyDict(
            window=window,
            column=Column(max_size=max_size, min_size=min_size),
            instance=None,
        )
        self.wx_object = self._state.column.wx_object
        self.refresh()
    
    def refresh(self):
        self._state.window.SetSizer(self._state.column.wx_object)
        self._state.window.Layout()
        self._state.window.Center(wx.BOTH)
    
    def add(self, *args, **kwargs):
        self._state.column.add(*args, **kwargs)
    
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        self.refresh()
        
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)