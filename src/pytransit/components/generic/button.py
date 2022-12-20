from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename
from pytransit.specific_tools import  gui_tools

class Button:
    """
        Overview:
            self.wx_object
            self.events.on_click
    """
    def __init__(self, text="", background_color=None, default_size=(250, -1), min_size=None, max_size=None):
        frame       = gui.frame
        
        # 
        # wx_object
        # 
        def wx_object(parent):
            if self._state.instance:
                return self._state.instance
            else:
                delayed_wx_object = GenBitmapTextButton(
                    frame,
                    -1,
                    gui_tools.bit_map,
                    text,
                    size=wx.Size(*default_size),
                )
                if background_color: delayed_wx_object.SetBackgroundColour(background_color)
                if max_size        : delayed_wx_object.SetMaxSize(wx.Size(*max_size))
                if min_size        : delayed_wx_object.SetMinSize(wx.Size(*min_size))
                
                delayed_wx_object.Bind(wx.EVT_BUTTON, lambda event: tuple(each(event) for each in self._state.callbacks))
            return delayed_wx_object
        
        self.wx_object = wx_object
        self.events = LazyDict(
            on_click=lambda func: self._state.callbacks.append(func),
        )
        self._state = LazyDict(
            instance=None,
            callbacks=[],
        )
        
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        self.Layout()
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)