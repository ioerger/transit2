from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools

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
        windowSizer = wx.BoxSizer(wx.VERTICAL)
        
        
        self.wx_object = wx_object
        self.events = LazyDict(
            # None
        )
        self._state = LazyDict(
            # None
        )
    
    def add(self, component, side, vertical_alignment, horizontal_alignment):
        # side:
            # wx.TOP
            # wx.BOTTOM
            # wx.LEFT
            # wx.RIGHT
            # wx.ALL

        idk_what_this_is = 1
        self.wx_object.Add(
            component.wx_object,
            idk_what_this_is,
            wx.ALL | wx.EXPAND,
            5
        )
        
        self.SetSizer(windowWrapper)
        self.Layout()
        
        self._state.index += 1
        wx_object.InsertItem(self._state.index, f"")
        if not isinstance(python_obj, dict):
            python_obj = python_obj.__dict__
        for each_key, each_value in python_obj.items():
            column_index = self._key_to_column_index(each_key)
            wx_object.SetItem(self._state.index, column_index, f"{each_value}")
    
    @property
    def length(self):
        return self._state.index+1
    
    @property
    def selected(self):
        selected = []
        current = -1
        while True:
            next = self.wx_object.GetNextSelected(current)
            if next == -1:
                break
            path = self.wx_object.GetItem(next)
            selected.append(path)
            current = next
        return selected
    
    # TODO: make a way to set the ones that are selected
    # @selected.setter
    # def selected(self, value):
    #     self._selected = value
    
    def __len__(self):
        return self.length
    
    def __enter__(self):
        return self.wx_object, self
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)