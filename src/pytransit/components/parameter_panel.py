import os

from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename, working_directory
import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
import pytransit.images as images
import pytransit

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table
                
# 
# options window
# 
panel = LazyDict()
def create_panel_area(_):
    panel.progress_count = 0
    
    # 
    # options window
    # 
    if True:
        panel.sizer = wx.BoxSizer(wx.VERTICAL)
        
        # 
        # Logo Section
        # 
        if True:
            logo_img = wx.StaticBitmap(
                universal.frame,
                wx.ID_ANY,
                wx.NullBitmap,
                wx.DefaultPosition,
                wx.DefaultSize,
                0,
            )
            panel.sizer.Add(logo_img, 0, wx.ALL | wx.ALIGN_CENTER, 5)
            logo_img.SetBitmap(images.transit_logo2.GetImage().ConvertToBitmap())
        
        # 
        # versionLabel
        # 
        if True:
            version_label = wx.StaticText(
                universal.frame,
                wx.ID_ANY,
                u"",
                wx.DefaultPosition,
                (100, 25),
                wx.ALIGN_CENTRE,
            )
            version_label.Wrap(-1)
            version_label.SetFont(wx.Font(10, 74, 90, 92, False, "Sans"))
            version_label.SetLabel(pytransit.__version__)
            

            panel.sizer.Add(
                version_label, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
            )
        
        # 
        # methodInfoSizer
        # 
        if True:
            panel.method_info_text = wx.StaticBox(universal.frame, wx.ID_ANY, u"Instructions")
            panel.method_info_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            panel.method_info_sizer = wx.StaticBoxSizer(panel.method_info_text, wx.VERTICAL)
            
            # 
            # methodShortText
            # 
            if True:
                panel.method_short_text = wx.StaticText(
                    universal.frame, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
                )
                panel.method_short_text.Wrap(250)
                panel.method_short_text.Hide()
                panel.method_info_sizer.Add(
                    panel.method_short_text, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                )
            
            # 
            # methodLongText
            # 
            if True:
                panel.method_long_text = wx.StaticText(
                    universal.frame, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
                )
                panel.method_long_text.Wrap(250)
                panel.method_long_text.Hide()
                panel.method_info_sizer.Add(
                    panel.method_long_text, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                )
            
            # 
            # methodDescText
            # 
            if True:

                panel.method_desc_text = wx.StaticText(
                    universal.frame, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
                )
                panel.method_desc_text.Wrap(250)
                panel.method_desc_text.Hide()
                panel.method_info_sizer.Add(
                    panel.method_desc_text, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                )
            
            # 
            # methodTnText
            # 
            if True:
                panel.method_tn_text = wx.StaticText(
                    universal.frame, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
                )
                panel.method_tn_text.Wrap(250)

                font = wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.BOLD)
                panel.method_tn_text.SetFont(font)
                panel.method_tn_text.Hide()

                panel.method_info_sizer.Add(
                    panel.method_tn_text, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                )
            
            # 
            # methodInstructions
            # 
            if True:
                panel.method_instructions = wx.StaticText(
                    universal.frame,
                    wx.ID_ANY,
                    universal.frame.instructions_text,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                panel.method_instructions.Wrap(250)
                panel.method_info_sizer.Add(
                    panel.method_instructions, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                )
            
            panel.sizer.Add(panel.method_info_sizer, 0, wx.ALL | wx.EXPAND, 5)
        
        # 
        # Method Options
        # 
        if True:
            panel.method_sizer_text = wx.StaticBox(universal.frame, wx.ID_ANY, u"Method Options")
            panel.method_sizer_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            panel.method_sizer = wx.StaticBoxSizer(panel.method_sizer_text, wx.VERTICAL)
            
            # 
            # methodPanel1
            # 
            if True:
                panel.method_panel1 = wx.Panel(
                    universal.frame,
                    wx.ID_ANY,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    wx.TAB_TRAVERSAL,
                )
                panel.method_panel1.SetMinSize(wx.Size(50, 1))
                panel.method_sizer.Add(panel.method_panel1, 0, wx.ALL, 5)
            
        panel.sizer.Add(panel.method_sizer, 0, wx.EXPAND, 5)

    
    # progress
    panel.progress_panel = wx.Panel(
        universal.frame,
        wx.ID_ANY,
        wx.DefaultPosition,
        wx.DefaultSize,
        wx.TAB_TRAVERSAL,
    )
    if True:
        progress_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # 
        # text
        # 
        if True:
            panel.progress_label = wx.StaticText(
                panel.progress_panel,
                wx.ID_ANY,
                u"Progress",
                wx.DefaultPosition,
                wx.DefaultSize,
                0,
            )
            panel.progress_label.Wrap(-1)
            progress_sizer.Add(panel.progress_label, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        if True:
            panel.progress = wx.Gauge(
                panel.progress_panel,
                wx.ID_ANY,
                20,
                wx.DefaultPosition,
                wx.DefaultSize,
                wx.GA_HORIZONTAL | wx.GA_SMOOTH,
            )
            progress_sizer.Add(panel.progress, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

    # panel.progress_panel.BackgroundColour = (0, 0, 250)
    panel.progress_panel.SetSizer(progress_sizer)
    panel.progress_panel.SetMaxSize(wx.Size(100, 100))
    panel.progress_panel.Layout()

    # progress_sizer.Fit( panel.progress_panel )
    panel.method_sizer.Add(
        panel.progress_panel,
        0,
        wx.ALL | wx.ALIGN_CENTER_HORIZONTAL,
        5,
    )

    # panel.method_sizer.Add( panel.global_label, 1, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
    # panel.method_sizer.Hide()
    panel.progress.SetRange(50)
    #########

    hide_progress_section()
    panel.method_sizer_text.Hide()
    
    return panel.sizer
    
def hide_all_options():
    from pytransit.analysis import methods
    
    hide_progress_section()
    for name in methods:
        methods[name].gui.Hide()

def hide_progress_section():
    panel.progress_label.Hide()
    panel.progress.Hide()

def show_progress_section():
    panel.progress_label.Show()
    panel.progress.Show()

def update_progress(msg):
    """"""
    method, count = msg
    panel.progress_count = count
    try:
        panel.progress.SetValue(panel.progress_count)
    except:
        pass

def set_progress_range(count):
    with gui_tools.nice_error_log:
        panel.progress.SetRange(count)

def finish_run():
    with gui_tools.nice_error_log:
        panel.progress_count = 0
        panel.progress.SetValue(panel.progress_count)