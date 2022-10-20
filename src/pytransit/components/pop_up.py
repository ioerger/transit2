from pytransit.specific_tools.transit_tools import wx
from pytransit.globals import gui, cli, root_folder, debugging_enabled

def create_pop_up(parent_wx_object):
    """
    @create_pop_up(parent_panel)
    def create_pop_up_contents(pop_up_panel, sizer, refresh, close):
        print("test")
    """
    def wrapper(function_being_wrapped):
        dialog_window = wx.Dialog(parent_wx_object, wx.FRAME_FLOAT_ON_PARENT)
        pop_up_sizer = wx.BoxSizer(wx.VERTICAL)
        dialog_window.SetSizer(pop_up_sizer)
        
        modal_status = None
        def refresh():
            nonlocal modal_status
            dialog_window.Layout()
            pop_up_sizer.Fit(dialog_window)
            modal_status = dialog_window.ShowModal()
        
        def close():
            nonlocal modal_status
            dialog_window.Destroy()
        
        function_being_wrapped(dialog_window, pop_up_sizer, refresh, close)

        dialog_window.Layout()
        pop_up_sizer.Fit(dialog_window)
        selected_path = dialog_window.ShowModal()

        return function_being_wrapped
    
    return wrapper