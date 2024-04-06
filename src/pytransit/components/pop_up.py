from pytransit.specific_tools.transit_tools import wx
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools import gui_tools

default_popup_position = (50, 50)

def create_pop_up(parent_wx_object, *, min_width=None, min_height=None, title=""):
    """
    @create_pop_up(parent_panel)
    def create_pop_up_contents(pop_up_panel, sizer, refresh, close):
        print("test")
    """
    def wrapper(function_being_wrapped):
        dialog_window = wx.Dialog(parent_wx_object, id=wx.ID_ANY, title=title, pos=default_popup_position, size=wx.Size((min_width or -1, min_height or -1)), style=wx.DEFAULT_DIALOG_STYLE | wx.MAXIMIZE_BOX | wx.STAY_ON_TOP | wx.RESIZE_BORDER, name="")
        
        pop_up_sizer = wx.BoxSizer(wx.VERTICAL)
        dialog_window.SetSizer(pop_up_sizer)
        
        modal_status = 1
        def refresh(*args):
            nonlocal modal_status
            pop_up_sizer.Fit(dialog_window)
            dialog_window.Layout()
            # set min size because pop_up_sizer.SetMinSize doesn't actually do its job
            fitted_width, fitted_height = dialog_window.GetSize()
            dialog_window.SetSize((
                max((min_width or -1), fitted_width),
                max((min_height or -1), fitted_height)
            ))
            dialog_window.Layout()
            if modal_status==1: modal_status = 2; dialog_window.ShowModal() 
        
        def close(*args):
            nonlocal modal_status
            modal_status = False
            dialog_window.Destroy()
        
        with gui_tools.nice_error_log:
            function_being_wrapped(dialog_window, pop_up_sizer, refresh, close)

        refresh()

        return function_being_wrapped
    
    return wrapper
