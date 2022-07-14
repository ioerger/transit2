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
frame = None
def create_panel_area(frame_input):
    
    global frame
    frame = frame_input
    
    
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
                frame,
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
                frame,
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
            panel.method_info_text = wx.StaticBox(frame, wx.ID_ANY, u"Instructions")
            panel.method_info_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            panel.method_info_sizer = wx.StaticBoxSizer(panel.method_info_text, wx.VERTICAL)
            
            # 
            # methodShortText
            # 
            if True:
                panel.method_short_text = wx.StaticText(
                    frame, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
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
                    frame, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
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
                    frame, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
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
                    frame, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
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
                    frame,
                    wx.ID_ANY,
                    frame.instructions_text,
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
            panel.method_sizer_text = wx.StaticBox(frame, wx.ID_ANY, u"Method Options")
            panel.method_sizer_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            panel.method_sizer = wx.StaticBoxSizer(panel.method_sizer_text, wx.VERTICAL)
            
            # 
            # methodPanel1
            # 
            if True:
                panel.method_panel1 = wx.Panel(
                    frame,
                    wx.ID_ANY,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    wx.TAB_TRAVERSAL,
                )
                panel.method_panel1.SetMinSize(wx.Size(50, 1))
                panel.method_sizer.Add(panel.method_panel1, 0, wx.ALL, 5)
            
            # 
            # globalLabel
            # 
            if True:
                panel.global_label = wx.StaticText(
                    frame,
                    wx.ID_ANY,
                    u"Global Options",
                    wx.DefaultPosition,
                    (130, 20),
                    0,
                )
                panel.global_label.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
                panel.method_sizer.Add(panel.global_label, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)
            
            # 
            # globalPanel
            # 
            if True:
                panel.global_panel = wx.Panel(
                    frame,
                    wx.ID_ANY,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    wx.TAB_TRAVERSAL,
                )
                panel.global_panel.SetMinSize(wx.Size(280, 150))
                panel.global_panel.SetMaxSize(wx.Size(-1, -1))

                global_sizer_vt = wx.BoxSizer(wx.VERTICAL)
                n_term_sizer    = wx.BoxSizer(wx.HORIZONTAL)
                c_term_sizer    = wx.BoxSizer(wx.HORIZONTAL)

                # N TERMINUS - GLOBAL
                panel.global_n_terminus_label = wx.StaticText(
                    panel.global_panel,
                    wx.ID_ANY,
                    u"Ignore N-Terminus %:",
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                panel.global_n_terminus_label.Wrap(-1)
                panel.global_n_terminus_text = wx.TextCtrl(panel.global_panel, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0)
                panel.global_n_terminus_icon = pytransit.analysis.base.InfoIcon(
                    panel.global_panel,
                    wx.ID_ANY,
                    tooltip="Ignores a fraction of the ORF, beginning at the N-terminal end. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",
                )
                n_term_sizer.Add(panel.global_n_terminus_label, 1, wx.ALIGN_CENTER, 5)
                n_term_sizer.Add(panel.global_n_terminus_text , 1, wx.ALIGN_CENTER, 5)
                n_term_sizer.Add(panel.global_n_terminus_icon , 1, wx.ALIGN_CENTER, 5)

                # C TERMINUS - GLOBAL
                panel.global_c_terminus_label = wx.StaticText(
                    panel.global_panel,
                    wx.ID_ANY,
                    u"Ignore C-Terminus %:",
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                panel.global_c_terminus_label.Wrap(-1)
                panel.global_c_terminus_text = wx.TextCtrl(panel.global_panel, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0)
                panel.global_c_terminus_icon = pytransit.analysis.base.InfoIcon(
                    panel.global_panel,
                    wx.ID_ANY,
                    tooltip="Ignores a fraction of the ORF, beginning at the C-terminal end. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",
                )

                c_term_sizer.Add(panel.global_c_terminus_label, 1, wx.ALIGN_CENTER_VERTICAL, 5)
                c_term_sizer.Add(panel.global_c_terminus_text, 1, wx.ALIGN_CENTER_VERTICAL, 5)
                c_term_sizer.Add(panel.global_c_terminus_icon, 1, wx.ALIGN_CENTER, 5)

                # Control Libraries text - GLOBAL
                ctrl_lib_sizer = wx.BoxSizer(wx.HORIZONTAL)
                panel.ctrl_lib_label = wx.StaticText(
                    panel.global_panel,
                    wx.ID_ANY,
                    u"Control Libraries:",
                    wx.DefaultPosition,
                    (170, -1),
                    0,
                )
                panel.ctrl_lib_label.Wrap(-1)
                panel.ctrl_lib_text = wx.TextCtrl(panel.global_panel, wx.ID_ANY, "", wx.DefaultPosition, (-1, -1), 0)
                panel.ctrl_lib_tip = pytransit.analysis.base.InfoIcon(
                    panel.global_panel,
                    wx.ID_ANY,
                    tooltip="String of letters representing an \
                    identifier for the libraries the datasets belong to. For example, if adding three \
                    datasets of different libraries, change the string to 'ABC'. Set of letters used  \
                    must match those in Experimental datasets. Keep empty or with all letters equal, e.g. \
                    'AAA', to do regular resampling.",
                )

                panel.ctrl_lib_text.Disable()
                ctrl_lib_sizer.Add(panel.ctrl_lib_label, 0, wx.ALIGN_CENTER_VERTICAL, 5)
                ctrl_lib_sizer.Add(panel.ctrl_lib_text, 0, wx.ALIGN_CENTER_VERTICAL, 5)
                ctrl_lib_sizer.Add(panel.ctrl_lib_tip, 0, wx.ALIGN_CENTER_VERTICAL, 5)

                # Experimental Libraries text - GLOBAL
                exp_lib_sizer = wx.BoxSizer(wx.HORIZONTAL)
                panel.exp_lib_label = wx.StaticText(
                    panel.global_panel,
                    wx.ID_ANY,
                    u"Experimental Libraries:",
                    wx.DefaultPosition,
                    (170, -1),
                    0,
                )
                panel.exp_lib_label.Wrap(-1)
                panel.exp_lib_text = wx.TextCtrl(
                    panel.global_panel, wx.ID_ANY, "", wx.DefaultPosition, (-1, -1), 0
                )
                panel.exp_lib_tip = pytransit.analysis.base.InfoIcon(
                    panel.global_panel,
                    wx.ID_ANY,
                    tooltip="String of letters representing an identifier for the libraries the datasets \
                    belong to. For example, if adding three datasets of different libraries, change the \
                    string to 'ABC'. Set  of letters used must match those in Control datasets. Keep \
                    empty or with all letters equal, e.g. 'AAA', to do regular resampling.",
                )

                panel.exp_lib_text.Disable()
                exp_lib_sizer.Add(panel.exp_lib_label, 0, wx.ALIGN_CENTER_VERTICAL, 5)
                exp_lib_sizer.Add(panel.exp_lib_text , 0, wx.ALIGN_CENTER_VERTICAL, 5)
                exp_lib_sizer.Add(panel.exp_lib_tip  , 0, wx.ALIGN_CENTER_VERTICAL, 5)

                global_sizer_vt.Add(n_term_sizer  , 1, wx.EXPAND, 5)
                global_sizer_vt.Add(c_term_sizer  , 1, wx.EXPAND, 5)
                global_sizer_vt.Add(ctrl_lib_sizer, 1, wx.EXPAND, 5)
                global_sizer_vt.Add(exp_lib_sizer , 1, wx.EXPAND, 5)

                panel.global_panel.SetSizer(global_sizer_vt)
                panel.global_panel.Layout()
                global_sizer_vt.Fit(panel.global_panel)
                panel.method_sizer.Add(panel.global_panel, 1, wx.ALIGN_CENTER_HORIZONTAL, 5)
            
        panel.sizer.Add(panel.method_sizer, 0, wx.EXPAND, 5)

    
    # progress
    panel.progress_panel = wx.Panel(
        frame,
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

    # options_window.Fit()

    hide_progress_section()
    hide_global_options()
    panel.method_sizer_text.Hide()
    
    return panel.sizer
    

def method_select_func(selected_name, event):
    # If empty is selected
    if selected_name == "[Choose Method]":
        method_wrap_width = 250
        hide_all_options()
        panel.method_info_text.SetLabel(u"Instructions")
        panel.method_instructions.Show()
        panel.method_instructions.SetLabel(frame.instructions_text)
        panel.method_instructions.Wrap(method_wrap_width)
        panel.method_short_text.Hide()
        panel.method_long_text.Hide()
        panel.method_tn_text.Hide()
        panel.method_desc_text.Hide()

        panel.method_choice = ""
    else:
        show_global_options()
        panel.method_sizer_text.Show()
        
        from pytransit.analysis import methods
        
        matched_name = None
        # Get selected Method and hide Others
        for name in methods:
            methods[name].gui.Hide()
            methods[name].gui.GlobalHide()
            methods[name].gui.GlobalDisable()

            if methods[name].fullname() == selected_name:
                matched_name = name
        
        print(f'''selected_name = {selected_name}''')
        print(f'''methods = {methods}''')
        if matched_name in methods:
            name = matched_name
            panel.method_info_text.SetLabel("%s" % methods[name].long_name)

            panel.method_tn_text.Show()
            panel.method_tn_text.SetLabel(methods[name].getTransposonsText())
            panel.method_tn_text.Wrap(250)

            panel.method_desc_text.Show()
            panel.method_desc_text.SetLabel(methods[name].getDescriptionText())
            panel.method_desc_text.Wrap(250)
            panel.method_instructions.SetLabel(" ")
            methods[name].gui.Show()
            methods[name].gui.Show()
            methods[name].gui.GlobalEnable()
            frame.statusBar.SetStatusText("[%s]" % methods[name].short_name)

        show_progress_section()
        panel.method_choice = selected_name

    frame.Layout()
    if frame.verbose:
        transit_tools.transit_message("Selected Method: %s" % (selected_name))


def hide_all_options():
    from pytransit.analysis import methods
    
    hide_global_options()
    hide_progress_section()
    for name in methods:
        methods[name].gui.Hide()

def hide_global_options():
    panel.global_label.Hide()
    hide_global_c_terminus()
    hide_global_n_terminus()
    hide_global_libraries()

def hide_global_c_terminus():
    panel.global_c_terminus_label.Hide()
    panel.global_c_terminus_text.Hide()
    panel.global_c_terminus_icon.Hide()

def hide_global_n_terminus():
    panel.global_n_terminus_label.Hide()
    panel.global_n_terminus_text.Hide()
    panel.global_n_terminus_icon.Hide()

def hide_global_libraries():
    hide_global_ctrl_libraries()   # BOOKMARK: control/experiment fragments
    hide_global_exp_libraries()

def hide_global_ctrl_libraries():
    panel.ctrl_lib_label.Hide()
    panel.ctrl_lib_text.Hide()
    panel.ctrl_lib_tip.Hide()

def hide_global_exp_libraries():
    panel.exp_lib_label.Hide()
    panel.exp_lib_text.Hide()
    panel.exp_lib_tip.Hide()

def show_global_options():
    panel.global_label.Show()
    show_global_c_terminus()
    show_global_n_terminus()
    show_global_libraries()

def show_global_c_terminus():
    panel.global_c_terminus_label.Show()
    panel.global_c_terminus_text.Show()
    panel.global_c_terminus_icon.Show()

def show_global_n_terminus():
    panel.global_n_terminus_text.Show()
    panel.global_n_terminus_label.Show()
    panel.global_n_terminus_icon.Show()

def show_global_libraries():
    show_global_ctrl_libraries()
    show_global_exp_libraries()

def show_global_ctrl_libraries():
    panel.ctrl_lib_label.Show()
    panel.ctrl_lib_text.Show()
    panel.ctrl_lib_tip.Show()

def show_global_exp_libraries():
    panel.exp_lib_label.Show()
    panel.exp_lib_text.Show()
    panel.exp_lib_tip.Show()

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