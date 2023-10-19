import os
import sys

from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.generic_tools import misc
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename
from pytransit.specific_tools import  gui_tools, transit_tools
import pytransit

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table

# 
# options window
# 
@misc.singleton
class panel:
    progress_sizer = None
    
    initial_instructions_text = """
        1. Choose the annotation file ("prot table") that corresponds to the datasets to be analyzed.
        
        2. Click "Load Combined Wig & Metadata", add both files
        
        3. (Optional) If you wish to visualize their read counts, select the desired datasets and click on the "View" button.
        
        4. Select the desired analysis method from the dropdown menu on the top-right of the window, and follow its instructions.
    """.replace("\n        ","\n")
    
    @property
    def max_width(self): # this is what the width SHOULD be, but wx does not always make it as such
        return int(gui.width * 0.3)
    
    @property
    def instructions_height(self): # this is what the width SHOULD be, but wx does not always make it as such
        return int(gui.height*0.21)
    
    parameter_proportion = 3
    parameter_padding = 15 # doesnt seem to do anything
    
    def refresh(self):
        gui.frame.Refresh()
        
def create_panel_area(_):
    import pytransit.components.images as images
    panel.progress_percent = 0
    
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
                gui.frame,
                wx.ID_ANY,
                wx.NullBitmap,
                wx.DefaultPosition,
                wx.DefaultSize,
                0,
            )
            panel.sizer.Add(logo_img, proportion=0, flag=wx.ALL | wx.ALIGN_CENTER, border=5)
            logo_img.SetBitmap(images.transit_logo2.GetImage().ConvertToBitmap())
        
        # 
        # versionLabel
        # 
        if True:
            version_label = wx.StaticText(
                gui.frame,
                wx.ID_ANY,
                "",
                wx.DefaultPosition,
                (100, 25),
                wx.ALIGN_CENTRE,
            )
            version_label.Wrap(panel.max_width)
            version_label.SetFont(wx.Font(10, 74, 90, 92, False, "Sans"))
            version_label.SetLabel(pytransit.__version__)
            

            panel.sizer.Add(
                version_label, proportion=0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=5
            )
        
        # 
        # methodInfoSizer
        # 
        if True:
            panel.method_info_text = wx.StaticBox(gui.frame, wx.ID_ANY, "Instructions")
            panel.method_info_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            panel.method_info_sizer = wx.StaticBoxSizer(panel.method_info_text, wx.VERTICAL)
            
            # 
            # Method Name
            # 
            if True:
                panel.method_name = wx.StaticText(
                    gui.frame, wx.ID_ANY, "  ", wx.DefaultPosition, wx.DefaultSize, 0
                )
                panel.method_name.Show()
                panel.method_name.Wrap(panel.max_width)
                panel.method_info_sizer.Add(
                    panel.method_name, proportion=0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=0,
                )
            
           
            # 
            # Instructions
            # 
            if True:
                panel.method_instructions = wx.TextCtrl(
                    gui.frame,
                    size= wx.Size(panel.max_width, panel.instructions_height),
                    style= wx.TE_MULTILINE | wx.TE_READONLY | wx.EXPAND,
                )
                panel.method_instructions.SetMinSize(wx.Size(panel.max_width, panel.instructions_height))
                panel.method_instructions.SetMaxSize(wx.Size(panel.max_width, panel.instructions_height))
                panel.method_instructions.SetValue(panel.initial_instructions_text)
                panel.method_info_sizer.Add(
                    panel.method_instructions, proportion=0, flag=wx.ALL | wx.EXPAND, border=5
                )
            
            panel.sizer.Add(panel.method_info_sizer, proportion=0, flag=wx.ALL | wx.EXPAND, border=0)
        
        # 
        # Method Options
        # 
        if True:
            panel.method_sizer = wx.BoxSizer(wx.VERTICAL)
            
        panel.sizer.Add(
            panel.method_sizer,
            proportion=int(panel.parameter_proportion),
            flag=wx.EXPAND,
            border=panel.parameter_padding,
        )
        
        # 
        # progress bar
        # 
        if True:
            progress_sizer = wx.BoxSizer(wx.VERTICAL)
            
            # 
            # text
            # 
            if True:
                panel.progress_label = wx.StaticText(
                    gui.frame,
                    wx.ID_ANY,
                    "Progress",
                    wx.DefaultPosition,
                    wx.Size(-1, -1),
                    0,
                )
                panel.progress_label.Wrap(panel.max_width)
                progress_sizer.Add(panel.progress_label, proportion=0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=5)

            if True:
                panel.progress = wx.Gauge(
                    gui.frame,
                    wx.ID_ANY,
                    20,
                    wx.DefaultPosition,
                    wx.Size(int(panel.max_width*0.7), 10),
                    wx.GA_HORIZONTAL | wx.GA_SMOOTH,
                )
                progress_sizer.Add(panel.progress, proportion=0, flag=wx.ALL | wx.EXPAND, border=0)

        panel.progress_sizer = progress_sizer
        panel.sizer.Add(panel.progress_sizer, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL, border=0)
        panel.sizer.Add(10,20) # vertical padding
        
    
    panel.progress.SetRange(1000)
    panel.progress_label.Hide()
    panel.progress.Hide()
    
    return panel.sizer

old_panel = None
def set_panel(new_panel):
    with gui_tools.nice_error_log:
        global old_panel
        if old_panel != None:
            old_panel.Hide()
        
        try: panel.method_sizer.Detach(new_panel)
        except Exception as error: print(error)
        
        panel.method_sizer.Add(
            new_panel,
            proportion=1,
            flag=wx.ALL|wx.EXPAND,
            border=gui_tools.default_padding,
        )
        new_panel.Show()
        
        panel.refresh()
        panel.method_sizer.Fit(new_panel)
        #new_panel.SetBackgroundColour(gui_tools.color.light_gray)
        old_panel = new_panel
        panel.progress_label.Show()
        panel.progress.Show()
        panel.refresh()

def set_instructions( title_text, sub_text,  method_specific_instructions,):
    with gui_tools.nice_error_log:
        panel.method_info_text.SetLabel(title_text +" Instructions:")
        panel.method_info_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        panel.method_info_text.Show()
        
        panel.method_name.SetLabel("--"+sub_text+"--")
        #panel.method_name.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        panel.method_name.Wrap(panel.max_width)
        
        
        panel.method_instructions.SetValue(method_specific_instructions)
        panel.method_instructions.Show()

def progress_update(text, percent):
    string = f" {text}   \r"
    
    # update current line
    sys.stdout.write(string)
    sys.stdout.flush()
    
    if HAS_WX:
        from pytransit.specific_tools import gui_tools
        # update progress bar
        panel.progress_percent = percent
        thousands = round(panel.progress_percent*10)
        try:
            panel.progress.SetValue(thousands)
        except:
            pass
        
        if gui.is_active:
            # update status text
            gui_tools.set_status(string)
            
            wx.Yield() # to get the UI to update