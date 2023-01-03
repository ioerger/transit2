from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename
from pytransit.specific_tools import  gui_tools
from pytransit.components.generic.box import Column, Row

class InnerFrame:
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, parent, title, default_size=(1350, 975), min_size=None, max_size=None):
        frame = gui.frame
        
        # 
        # wx_object
        # 
        wx.Frame.__init__(
            frame,
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
            frame=frame,
            column=Column(max_size=max_size, min_size=min_size),
            instance=None,
        )
        self.wx_object = self._state.column.wx_object
        self.refresh()
    
    def refresh(self):
        self._state.frame.SetSizer(self._state.column.wx_object)
        self._state.frame.Layout()
        self._state.frame.Center(wx.BOTH)
    
    def add(self, *args, **kwargs):
        self._state.column.add(*args, **kwargs)
    
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        self.refresh()
        
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)