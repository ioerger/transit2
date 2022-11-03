from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename
from pytransit.specific_tools import logging, gui_tools

class Table:
    """
        Overview:
            self.wx_object
            self.events.on_select   # decorator (use @)
            self.selected           # list of wx objects, TODO: have it return the python_obj
            self.length
            self.add(python_obj)
    """
    def __init__(self, initial_columns=None, column_width=None, min_size=(1,1), max_size=(-1, -1), frame=None):
        frame        = gui.frame if not frame else frame
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
        self.min_size = min_size
        self.max_size = max_size
        self.wx_object = wx_object
        
        self.refresh_size()
        self.wx_object.InsertColumn(0, "", width=0) # first one is some kind of special name. Were going to ignore it    
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
    
    def refresh_size(self):
        min_width, min_height = self.min_size
        max_width, max_height = self.max_size
        # set min size because pop_up_sizer.SetMinSize doesn't actually do its job
        fitted_width, fitted_height = self.wx_object.GetSize()
        if fitted_width  > max_width : fitted_width  = max_width
        if fitted_height > max_height: fitted_height = max_height
        if fitted_width  < min_width : fitted_width  = min_width
        if fitted_height < min_height: fitted_height = min_height
        self.wx_object.SetSize((
            int(fitted_width),
            int(fitted_height),
        ))
        self.wx_object.SetMinSize((
            int(fitted_width),
            int(fitted_height),
        ))
        self.wx_object.SetMaxSize((
            int(fitted_width),
            int(fitted_height),
        ))
    
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
            if not isinstance(each_key, str):
                continue
            if each_key.startswith("__"):
                continue
            column_index = self._key_to_column_index(each_key)
            self.wx_object.SetItem(self._state.index, column_index, f"{each_value}")
        
        self._state.data_values.append(python_obj)
    
    @property
    def length(self):
        return self._state.index+1
    
    @property
    def selected_wx_objects(self):
        selected = []
        current_selected_index = -1
        while True:
            next_selected_index = self.wx_object.GetNextSelected(current_selected_index)
            if next_selected_index == -1:
                break
            selected.append(
                self.wx_object.GetItem(next_selected_index)
            )
            current_selected_index = next_selected_index
        return selected
    
    @property
    def selected_rows(self):
        selected = []
        current_selected_index = -1
        while True:
            next_selected_index = self.wx_object.GetNextSelected(current_selected_index)
            if next_selected_index == -1:
                break
            selected.append(self.rows[next_selected_index])
            current_selected_index = next_selected_index
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