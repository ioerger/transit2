import os

from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename, working_directory
from pytransit.analysis import methods
import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
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
    
    options_window = wx.ScrolledWindow(
        frame,
        wx.ID_ANY,
        wx.DefaultPosition,
        wx.Size(-1, -1),
        wx.HSCROLL | wx.VSCROLL | wx.EXPAND,
    )
    options_window.SetScrollRate(5, 5)
    options_window.SetMinSize(wx.Size(310, 1000))
    
    # 
    # box
    # 
    if True:
        options_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # 
        # Logo Section
        # 
        if True:
            logo_img = wx.StaticBitmap(
                options_window,
                wx.ID_ANY,
                wx.NullBitmap,
                wx.DefaultPosition,
                wx.DefaultSize,
                0,
            )
            options_sizer.Add(logo_img, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        
        # 
        # versionLabel
        # 
        if True:
            version_label = wx.StaticText(
                options_window,
                wx.ID_ANY,
                u"",
                wx.DefaultPosition,
                (100, 25),
                wx.ALIGN_CENTRE,
            )
            version_label.Wrap(-1)
            version_label.SetFont(wx.Font(10, 74, 90, 92, False, "Sans"))

            options_sizer.Add(
                version_label, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
            )
        
        # 
        # methodInfoSizer
        # 
        if True:
            panel.method_info_text = wx.StaticBox(options_window, wx.ID_ANY, u"Instructions")
            panel.method_info_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            panel.method_info_sizer = wx.StaticBoxSizer(panel.method_info_text, wx.VERTICAL)
            
            # 
            # methodShortText
            # 
            if True:
                panel.method_short_text = wx.StaticText(
                    options_window, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
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
                    options_window, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
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
                    options_window, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
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
                    options_window, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0
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
                    options_window,
                    wx.ID_ANY,
                    self.instructions_text,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                panel.method_instructions.Wrap(250)
                panel.method_info_sizer.Add(
                    panel.method_instructions, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
                )
            
            options_sizer.Add(panel.method_info_sizer, 0, wx.ALL | wx.EXPAND, 5)
        
        # 
        # Method Options
        # 
        if True:
            panel.method_sizer_text = wx.StaticBox(options_window, wx.ID_ANY, u"Method Options")
            panel.method_sizer_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            panel.method_sizer = wx.StaticBoxSizer(panel.method_sizer_text, wx.VERTICAL)
            
            # 
            # methodPanel1
            # 
            if True:
                panel.method_panel1 = wx.Panel(
                    options_window,
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
                    options_window,
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
                    options_window,
                    wx.ID_ANY,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    wx.TAB_TRAVERSAL,
                )
                panel.global_panel.SetMinSize(wx.Size(280, 150))
                panel.global_panel.SetMaxSize(wx.Size(-1, -1))

                globalSizerVT = wx.BoxSizer(wx.VERTICAL)
                nTermSizer    = wx.BoxSizer(wx.HORIZONTAL)
                cTermSizer    = wx.BoxSizer(wx.HORIZONTAL)

                # N TERMINUS - GLOBAL
                self.globalNTerminusLabel = wx.StaticText(
                    panel.global_panel,
                    wx.ID_ANY,
                    u"Ignore N-Terminus %:",
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                self.globalNTerminusLabel.Wrap(-1)
                self.globalNTerminusText = wx.TextCtrl(panel.global_panel, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0)
                self.globalNTerminusIcon = pytransit.analysis.base.InfoIcon(
                    panel.global_panel,
                    wx.ID_ANY,
                    tooltip="Ignores a fraction of the ORF, beginning at the N-terminal end. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",
                )
                nTermSizer.Add(self.globalNTerminusLabel, 1, wx.ALIGN_CENTER, 5)
                nTermSizer.Add(self.globalNTerminusText , 1, wx.ALIGN_CENTER, 5)
                nTermSizer.Add(self.globalNTerminusIcon , 1, wx.ALIGN_CENTER, 5)

                # C TERMINUS - GLOBAL
                self.globalCTerminusLabel = wx.StaticText(
                    panel.global_panel,
                    wx.ID_ANY,
                    u"Ignore C-Terminus %:",
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                self.globalCTerminusLabel.Wrap(-1)
                self.globalCTerminusText = wx.TextCtrl(panel.global_panel, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0)
                self.globalCTerminusIcon = pytransit.analysis.base.InfoIcon(
                    panel.global_panel,
                    wx.ID_ANY,
                    tooltip="Ignores a fraction of the ORF, beginning at the C-terminal end. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",
                )

                cTermSizer.Add(self.globalCTerminusLabel, 1, wx.ALIGN_CENTER_VERTICAL, 5)
                cTermSizer.Add(self.globalCTerminusText, 1, wx.ALIGN_CENTER_VERTICAL, 5)
                cTermSizer.Add(self.globalCTerminusIcon, 1, wx.ALIGN_CENTER, 5)

                # Control Libraries text - GLOBAL
                ctrlLibSizer = wx.BoxSizer(wx.HORIZONTAL)
                self.ctrlLibLabel = wx.StaticText(
                    panel.global_panel,
                    wx.ID_ANY,
                    u"Control Libraries:",
                    wx.DefaultPosition,
                    (170, -1),
                    0,
                )
                self.ctrlLibLabel.Wrap(-1)
                self.ctrlLibText = wx.TextCtrl(panel.global_panel, wx.ID_ANY, "", wx.DefaultPosition, (-1, -1), 0)
                self.ctrlLibTip = pytransit.analysis.base.InfoIcon(
                    panel.global_panel,
                    wx.ID_ANY,
                    tooltip="String of letters representing an \
                    identifier for the libraries the datasets belong to. For example, if adding three \
                    datasets of different libraries, change the string to 'ABC'. Set of letters used  \
                    must match those in Experimental datasets. Keep empty or with all letters equal, e.g. \
                    'AAA', to do regular resampling.",
                )

                self.ctrlLibText.Disable()
                ctrlLibSizer.Add(self.ctrlLibLabel, 0, wx.ALIGN_CENTER_VERTICAL, 5)
                ctrlLibSizer.Add(self.ctrlLibText, 0, wx.ALIGN_CENTER_VERTICAL, 5)
                ctrlLibSizer.Add(self.ctrlLibTip, 0, wx.ALIGN_CENTER_VERTICAL, 5)

                # Experimental Libraries text - GLOBAL
                expLibSizer = wx.BoxSizer(wx.HORIZONTAL)
                self.expLibLabel = wx.StaticText(
                    panel.global_panel,
                    wx.ID_ANY,
                    u"Experimental Libraries:",
                    wx.DefaultPosition,
                    (170, -1),
                    0,
                )
                self.expLibLabel.Wrap(-1)
                self.expLibText = wx.TextCtrl(
                    panel.global_panel, wx.ID_ANY, "", wx.DefaultPosition, (-1, -1), 0
                )
                self.expLibTip = pytransit.analysis.base.InfoIcon(
                    panel.global_panel,
                    wx.ID_ANY,
                    tooltip="String of letters representing an identifier for the libraries the datasets \
                    belong to. For example, if adding three datasets of different libraries, change the \
                    string to 'ABC'. Set  of letters used must match those in Control datasets. Keep \
                    empty or with all letters equal, e.g. 'AAA', to do regular resampling.",
                )

                self.expLibText.Disable()
                expLibSizer.Add(self.expLibLabel, 0, wx.ALIGN_CENTER_VERTICAL, 5)
                expLibSizer.Add(self.expLibText , 0, wx.ALIGN_CENTER_VERTICAL, 5)
                expLibSizer.Add(self.expLibTip  , 0, wx.ALIGN_CENTER_VERTICAL, 5)

                globalSizerVT.Add(nTermSizer  , 1, wx.EXPAND, 5)
                globalSizerVT.Add(cTermSizer  , 1, wx.EXPAND, 5)
                globalSizerVT.Add(ctrlLibSizer, 1, wx.EXPAND, 5)
                globalSizerVT.Add(expLibSizer , 1, wx.EXPAND, 5)

                panel.global_panel.SetSizer(globalSizerVT)
                panel.global_panel.Layout()
                globalSizerVT.Fit(panel.global_panel)
                panel.method_sizer.Add(panel.global_panel, 1, wx.ALIGN_CENTER_HORIZONTAL, 5)
            
        options_sizer.Add(panel.method_sizer, 0, wx.EXPAND, 5)

    options_window.SetSizer(options_sizer)
    options_window.Layout()

    options_window.Fit()



def MethodSelectFunc(selected_name):
    # If empty is selected
    if selected_name == "[Choose Method]":
        method_wrap_width = 250
        HideAllOptions()
        panel.method_info_text.SetLabel(u"Instructions")
        panel.method_instructions.Show()
        panel.method_instructions.SetLabel(self.instructions_text)
        panel.method_instructions.Wrap(method_wrap_width)
        panel.method_short_text.Hide()
        panel.method_long_text.Hide()
        panel.method_tn_text.Hide()
        panel.method_desc_text.Hide()

        panel.method_choice = ""
    else:
        ShowGlobalOptions()
        panel.method_sizer_text.Show()

        # Get selected Method and hide Others
        for name in methods:
            methods[name].gui.Hide()
            methods[name].gui.GlobalHide()
            methods[name].gui.GlobalDisable()

            if methods[name].fullname() == selected_name:
                matched_name = name

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

        ShowProgressSection()
        panel.method_choice = selected_name

    frame.Layout()
    if frame.verbose:
        transit_tools.transit_message("Selected Method: %s" % (selected_name))



def HideAllOptions(self):
    HideGlobalOptions()
    HideProgressSection()
    for name in methods:
        methods[name].gui.Hide()

def HideGlobalOptions(self):
    panel.global_label.Hide()
    self.HideGlobalCTerminus()
    self.HideGlobalNTerminus()
    self.HideGlobalLibraries()

def HideGlobalCTerminus(self):
    self.globalCTerminusLabel.Hide()
    self.globalCTerminusText.Hide()
    self.globalCTerminusIcon.Hide()

def HideGlobalNTerminus(self):
    self.globalNTerminusLabel.Hide()
    self.globalNTerminusText.Hide()
    self.globalNTerminusIcon.Hide()

def HideGlobalLibraries(self):
    self.HideGlobalCtrlLibraries()
    self.HideGlobalExpLibraries()

def HideGlobalCtrlLibraries(self):
    self.ctrlLibLabel.Hide()
    self.ctrlLibText.Hide()
    self.ctrlLibTip.Hide()

def HideGlobalExpLibraries(self):
    self.expLibLabel.Hide()
    self.expLibText.Hide()
    self.expLibTip.Hide()

def ShowGlobalOptions(self):
    panel.global_label.Show()
    self.ShowGlobalCTerminus()
    self.ShowGlobalNTerminus()
    self.ShowGlobalLibraries()

def ShowGlobalCTerminus(self):
    self.globalCTerminusLabel.Show()
    self.globalCTerminusText.Show()
    self.globalCTerminusIcon.Show()

def ShowGlobalNTerminus(self):
    self.globalNTerminusText.Show()
    self.globalNTerminusLabel.Show()
    self.globalNTerminusIcon.Show()

def ShowGlobalLibraries(self):
    self.ShowGlobalCtrlLibraries()
    self.ShowGlobalExpLibraries()

def ShowGlobalCtrlLibraries(self):
    self.ctrlLibLabel.Show()
    self.ctrlLibText.Show()
    self.ctrlLibTip.Show()

def ShowGlobalExpLibraries(self):
    self.expLibLabel.Show()
    self.expLibText.Show()
    self.expLibTip.Show()

def HideProgressSection(self):
    self.progressLabel.Hide()
    self.progress.Hide()

def ShowProgressSection(self):
    self.progressLabel.Show()
    self.progress.Show()