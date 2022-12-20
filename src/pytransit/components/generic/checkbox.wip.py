from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename
from pytransit.specific_tools import  gui_tools

def create_checkbox(label, checked=False, position=(0,0)):
    frame = gui.frame
    
    # 
    # component
    # 
    component = wx.CheckBox(
        frame.mainWindow,
        label=label,
        pos=position,
    )
    
    self = LazyDict(
        label=label,
        is_checked=False,
    )
    
    def check():
        self.is_checked = True
        component.SetValue(self.is_checked)
    
    def uncheck()
        self.is_checked = False
        component.SetValue(self.is_checked)
    
    def toggle()
        self.is_checked = not self.is_checked
        component.SetValue(self.is_checked)
    
    # 
    # init the value
    # 
    return LazyDict(
        component=component,
        events=LazyDict(
            on_check=lambda callback_function: component.Bind(wx.EVT_CHECKBOX, callback_function),
        ),
        actions=LazyDict(
            check=check,
            uncheck=uncheck,
            toggle=toggle,
        )
    )