from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools

def create_table(column_width=100):
    window = gui_tools.window
    
    # 
    # component
    # 
    component = wx.ListCtrl(
        window.mainWindow,
        wx.ID_ANY,
        wx.DefaultPosition,
        wx.DefaultSize,
        wx.LC_REPORT | wx.SUNKEN_BORDER,
    )
    component.SetMaxSize(wx.Size(-1, 200))
    
    self = LazyDict(
        key_to_column_index={},
        index=-1,
        column_width=column_width,
    )
    
    # 
    # init the columns
    # 
    component.InsertColumn(0, "", width=0) # first one is some kind of special name. Were going to ignore it
    def key_to_column_index(key):
        if key not in self.key_to_column_index:
            index = len(self.key_to_column_index)+1
            self.key_to_column_index[key] = index
            component.InsertColumn(index, key, width=self.column_width)
            return index
        else:
            return self.key_to_column_index[key]
            
    
    def add(element):
        self.index += 1
        component.InsertItem(self.index, f"")
        if not isinstance(element, dict):
            element = element.__dict__
        for each_key, each_value in element.items():
            column_index = key_to_column_index(each_key)
            component.SetItem(self.index, column_index, f"{each_value}")
    
    return LazyDict(
        component=component,
        events=LazyDict(
            on_select=gui_tools.bind_to(component, wx.EVT_LIST_ITEM_SELECTED),
        ),
        actions=LazyDict(
            add=add,
        )
    )