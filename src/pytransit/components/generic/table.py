from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools

class Table:
    """
        Overview:
            self.wx_object
            self.events.on_select   # decorator (use @)
            self.selected           # list of wx objects, TODO: have it return the python_obj
            self.length
            self.add(python_obj)
    """
    def __init__(self, initial_columns=None, column_width=None, max_size=(-1, 200)):
        frame        = universal.frame
        column_width = column_width if column_width is not None else 100
        
        # 
        # wx_object
        # 
        wx_object = wx.ListCtrl(
            frame,
            wx.ID_ANY,
            wx.DefaultPosition,
            wx.DefaultSize,
            wx.LC_REPORT | wx.SUNKEN_BORDER,
        )
        wx_object.SetMaxSize(wx.Size(*max_size))
        wx_object.InsertColumn(0, "", width=0) # first one is some kind of special name. Were going to ignore it
        
        self.wx_object = wx_object
        self.events = LazyDict(
            on_select=lambda func: wx_object.Bind(wx.EVT_LIST_ITEM_SELECTED, func),
        )
        
        self._state = LazyDict(
            index               = -1,
            key_to_column_index = {},
            column_width        = column_width,
            initial_columns     = initial_columns or [],
            data_values         = [],
        )
        
        # create the inital columns
        for each_key in self._state.initial_columns:
            self._key_to_column_index(each_key)
        
    
    def _key_to_column_index(self, key):
        if key not in self._state.key_to_column_index:
            index = len(self._state.key_to_column_index)+1
            self._state.key_to_column_index[key] = index
            self.wx_object.InsertColumn(index, key, width=self._state.column_width)
            return index
        else:
            return self._state.key_to_column_index[key]
    
    def add(self, python_obj):
        self._state.index += 1
        self.wx_object.InsertItem(self._state.index, f"")
        if not isinstance(python_obj, dict):
            python_obj = python_obj.__dict__
        for each_key, each_value in python_obj.items():
            column_index = self._key_to_column_index(each_key)
            self.wx_object.SetItem(self._state.index, column_index, f"{each_value}")
        
        self._state.data_values.append(dict(python_obj))
    
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
    
    @property
    def rows(self):
        return list(self._state.data_values)
    
    @property
    def column_names(self):
        return list(self._state.key_to_column_index.keys())
    
    # TODO: make a way to set the ones that are selected
    # @selected.setter
    # def selected(self, value):
    #     self._selected = value
    
    def __len__(self):
        return self.length
    
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)