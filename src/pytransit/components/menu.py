from collections import defaultdict
from functools import partial

from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, stat_tools
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
import pytransit

selected_export_menu_item = None
documentation_url = "http://saclab.tamu.edu/essentiality/transit/transit.html"

def create_menu(frame):
    # must imported inside the function to avoid circular import
    import pytransit.components.parameter_panel as parameter_panel
    from pytransit.specific_tools.transit_tools import wx
    
    global selected_export_menu_item
    menu_bar = wx.MenuBar(0)
    frame = gui.frame
    
    @gui.add_menu("File", "Exit")
    def on_exit(event):
        gui.frame.Close()
    
    # 
    # Automated menu creation
    # 
    if True:
        def recursive_create_sub_menu(remaining_hierarchy):
            parent_menu = wx.Menu()
            for each_name, each_sub_value in remaining_hierarchy.items():
                # base-case option
                if callable(each_sub_value):
                    on_click_function = each_sub_value
                    menu_item = wx.MenuItem(
                        parent_menu,
                        wx.ID_ANY,
                        f"&{each_name}",
                        wx.EmptyString,
                        wx.ITEM_NORMAL,
                    )
                    parent_menu.Append(menu_item)
                    frame.Bind(wx.EVT_MENU, on_click_function, id=menu_item.GetId())
                # more heirarchy
                elif isinstance(each_sub_value, dict):
                    parent_menu.AppendSubMenu(
                        recursive_create_sub_menu(each_sub_value),
                        f"&{each_name}"
                    )
            
            return parent_menu
        
        # top level menu items are an edgecase
        for each_name, each_sub_value in gui.menu_heirarchy.items():
            menu_bar.Append(
                recursive_create_sub_menu(each_sub_value),
                f"&{each_name}"
            )
    # 
    # Help Menu
    # 
    if True:
        help_menu = wx.Menu()
        
        # 
        # Documentation
        # 
        if True:
            documentation_option = wx.MenuItem(
                help_menu,
                wx.ID_ANY,
                "&Documentation",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            help_menu.Append(documentation_option)
            def when_documentation_clicked(event):
                with gui_tools.nice_error_log:
                    from pytransit.generic_tools.misc import open_url
                    open_url(documentation_url)
                    
            frame.Bind(wx.EVT_MENU, when_documentation_clicked, id=documentation_option.GetId())
        
        # 
        # About
        #
        if True: 
            about_option = wx.MenuItem(
                help_menu, wx.ID_ANY, "&About", wx.EmptyString, wx.ITEM_NORMAL
            )
            help_menu.Append(about_option)
            def when_about_option_clicked(event):
                with gui_tools.nice_error_log:
                    import pytransit.components.images as images
                    
                    description = """TRANSIT is a tool for analysing TnSeq data. It provides an easy to use graphical interface and access to several different analysis methods that allow the user to determine essentiality within a single condition as well as between two conditions.


                        If you need to cite this tool, please use the following reference:

                        DeJesus, M.A., Ambadipudi, C., Baker, R., Sassetti, C., and Ioerger, T.R. (2015). TRANSIT - a Software Tool for Himar1 TnSeq Analysis. PLOS Computational Biology, 11(10):e1004401


                    """.replace("\n                    ","\n")

                    licence = """
                        TRANSIT is free software: you can redistribute it and/or modify
                        it under the terms of the GNU General Public License as published by
                        the Free Software Foundation, either version 3 of the License.


                        TRANSIT is distributed in the hope that it will be useful,
                        but WITHOUT ANY WARRANTY; without even the implied warranty of
                        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
                        GNU General Public License for more details.

                        You should have received a copy of the GNU General Public License
                        along with TRANSIT.  If not, see <http://www.gnu.org/licenses/>.
                    """.replace("\n            ", "\n")

                    info = wx.adv.AboutDialogInfo()
                    info.SetIcon(images.transit_logo2.GetIcon())
                    # images.transit_logo2.GetImage().ConvertToBitmap()
                    info.SetName("TRANSIT")
                    info.SetVersion(pytransit.__version__)
                    info.SetDescription(description)
                    info.SetCopyright("(C) 2015\n Michael A. DeJesus\nThomas R. Ioerger")
                    info.SetWebSite("http://saclab.tamu.edu/essentiality/transit/")
                    info.SetLicence(licence)
                    info.AddDeveloper("Michael A. DeJesus")
                    info.AddDeveloper("Thomas R. Ioerger")
                    info.AddDeveloper("Chaitra Ambadipudi")
                    info.AddDeveloper("Richard Baker")
                    info.AddDeveloper("Christopher Sassetti")
                    info.AddDeveloper("Eric Nelson")
                    wx.adv.AboutBox(info)
            frame.Bind(wx.EVT_MENU, when_about_option_clicked, id=about_option.GetId())
            
            menu_bar.Append(help_menu, "&Help")
    
    frame.SetMenuBar(menu_bar)