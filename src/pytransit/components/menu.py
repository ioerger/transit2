import os
from functools import partial

from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename, working_directory
from pytransit.analysis   import methods
from pytransit.export     import methods as export_methods
from pytransit.convert    import methods as convert_methods
from pytransit.norm_tools import methods as norm_methods
import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
import pytransit.images as images
import pytransit

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table
                
                

# TODO: cleanup formatting of nested items

def create_menu(self):
    if True:
        menu_bar = wx.MenuBar(0)
        file_menu_item = wx.Menu()
        export_menu_item = wx.Menu()
        selected_export_menu_item = wx.Menu()

        # Selected datasets
        export_menu_item.AppendSubMenu(
            selected_export_menu_item, "Selected Datasets"
        )

        file_menu_item.AppendSubMenu(export_menu_item, "Export")

        self.convertMenuItem = wx.Menu()
        self.annotationConvertPTToPTTMenu = wx.MenuItem(
            self.convertMenuItem,
            wx.ID_ANY,
            "prot_table to PTT",
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )
        self.convertMenuItem.Append(self.annotationConvertPTToPTTMenu)

        self.annotationConvertPTToGFF3Menu = wx.MenuItem(
            self.convertMenuItem,
            wx.ID_ANY,
            "prot_table to GFF3",
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )
        self.convertMenuItem.Append(self.annotationConvertPTToGFF3Menu)

        self.annotationConvertPTTToPT = wx.MenuItem(
            self.convertMenuItem,
            wx.ID_ANY,
            "PTT to prot_table",
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )

        self.convertMenuItem.Append(self.annotationConvertPTTToPT)

        # self.annotationConvertGFF3ToPT = wx.MenuItem( self.convertMenuItem, wx.ID_ANY, "GFF3 to prot_table", wx.EmptyString, wx.ITEM_NORMAL )
        # self.convertMenuItem.Append( self.annotationConvertGFF3ToPT )
        file_menu_item.AppendSubMenu(self.convertMenuItem, "Convert")

        self.fileExitMenuItem = wx.MenuItem(
            file_menu_item, wx.ID_ANY, "&Exit", wx.EmptyString, wx.ITEM_NORMAL
        )
        file_menu_item.Append(self.fileExitMenuItem)
        menu_bar.Append(file_menu_item, "&File")

        self.viewMenuItem = wx.Menu()
        self.scatterMenuItem = wx.MenuItem(
            self.viewMenuItem,
            wx.ID_ANY,
            "&Scatter Plot",
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )

        self.viewMenuItem.Append(self.scatterMenuItem)

        self.trackMenuItem = wx.MenuItem(
            self.viewMenuItem, wx.ID_ANY, "&Track View", wx.EmptyString, wx.ITEM_NORMAL
        )

        self.viewMenuItem.Append(self.trackMenuItem)
        menu_bar.Append(self.viewMenuItem, "&View")

        #
        self.methodsMenuItem = wx.Menu()
        self.himar1MenuItem = wx.Menu()
        self.tn5MenuItem = wx.Menu()

        self.methodsMenuItem.AppendSubMenu(self.himar1MenuItem, "&Himar1 Methods")
        self.methodsMenuItem.AppendSubMenu(self.tn5MenuItem, "&Tn5 Methods")
        menu_bar.Append(self.methodsMenuItem, "&Analysis")

        self.SetMenuBar(menu_bar)

        self.helpMenuItem = wx.Menu()
        self.documentationMenuItem = wx.MenuItem(
            self.helpMenuItem,
            wx.ID_ANY,
            "&Documentation",
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )
        self.helpMenuItem.Append(self.documentationMenuItem)
        self.aboutMenuItem = wx.MenuItem(
            self.helpMenuItem, wx.ID_ANY, "&About", wx.EmptyString, wx.ITEM_NORMAL
        )
        self.helpMenuItem.Append(self.aboutMenuItem)

        menu_bar.Append(self.helpMenuItem, "&Help")

        self.statusBar = self.CreateStatusBar(1, wx.STB_SIZEGRIP, wx.ID_ANY)
    
    # 
    # Export Menu Items
    # 
    for name in export_methods:
        export_methods[name].gui.defineMenuItem(self, export_methods[name].label)
        tempMenuItem = export_methods[name].gui.menuitem
        self.selectedExportMenuItem.Append(tempMenuItem)

        self.Bind(
            wx.EVT_MENU,
            partial(self.ExportSelectFunc, export_methods[name].label),
            tempMenuItem,
        )

    # Convert Menu Items
    for name in convert_methods:
        convert_methods[name].gui.defineMenuItem(self, convert_methods[name].label)
        tempMenuItem = convert_methods[name].gui.menuitem
        self.convertMenuItem.Append(tempMenuItem)

        self.Bind(
            wx.EVT_MENU,
            partial(self.ConvertSelectFunc, convert_methods[name].label),
            tempMenuItem,
        )
