from collections import defaultdict
from functools import partial

from pytransit.tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, stat_tools
import pytransit.components.qc_display as qc_display
from pytransit.globals import gui, cli, root_folder, debugging_enabled
import pytransit

selected_export_menu_item = None
documentation_url = "http://saclab.tamu.edu/essentiality/transit/transit.html"

def create_menu(frame):
    # must imported inside the function to avoid circular import
    import pytransit.components.parameter_panel as parameter_panel
    from pytransit.tools.transit_tools import wx
    
    global selected_export_menu_item
    menu_bar = wx.MenuBar(0)
    frame = gui.frame
    
    # 
    # File Menu
    # 
    if True:
        file_menu = wx.Menu()
    
        # 
        # Exit
        # 
        if True:
            exit_option = wx.MenuItem( file_menu, wx.ID_ANY, "&Exit", wx.EmptyString, wx.ITEM_NORMAL )
            file_menu.Append(exit_option)
            def when_exit_clicked(event):
                if frame.verbose: logging.log("Exiting Transit")
                frame.Close()
            frame.Bind(wx.EVT_MENU, when_exit_clicked, id=exit_option.GetId())
        
        menu_bar.Append(file_menu, "&File")
    
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
                    from pytransit.basics.misc import open_url
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


# UNUSED 
def annotation_gff3_to_pt(event):
    with gui_tools.nice_error_log:
        annotation_path = gui.annotation_path
        default_file = transit_tools.fetch_name(annotation_path) + ".prot_table"
        # default_dir = os.path.dirname(os.path.realpath(__file__))
        default_dir = os.getcwd()

        if not annotation_path:
            # NOTE: was a popup
            logging.error("Error: No annotation file selected.")
        else:
            output_path = frame.SaveFile(default_dir, default_file)
            if not output_path:
                return
            if frame.verbose:
                logging.log(
                    "Converting annotation file from GFF3 format to prot_table format"
                )

            output = open(output_path, "w")
            with open(annotation_path) as file:
                for line in file:
                    if line.startswith("#"):
                        continue
                    tmp = line.strip().split("\t")
                    chr = tmp[0]
                    type = tmp[2]
                    start = int(tmp[3])
                    end = int(tmp[4])
                    length = ((end - start + 1) / 3) - 1
                    strand = tmp[6]
                    features = dict([tuple(f.split("=")) for f in tmp[8].split(";")])
                    if "ID" not in features:
                        continue
                    orf = features["ID"]
                    name = features.get("Name", "-")
                    if name == "-":
                        name = features.get("name", "-")

                    desc = features.get("Description", "-")
                    if desc == "-":
                        desc = features.get("description", "-")
                    if desc == "-":
                        desc = features.get("Desc", "-")
                    if desc == "-":
                        desc = features.get("desc", "-")

                    someID = "-"
                    someID2 = "-"
                    COG = "-"
                    output.write(
                        "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n"
                        % (
                            desc,
                            start,
                            end,
                            strand,
                            length,
                            someID,
                            someID2,
                            name,
                            orf,
                            COG,
                        )
                    )
            output.close()
            if frame.verbose:
                logging.log("Finished conversion")