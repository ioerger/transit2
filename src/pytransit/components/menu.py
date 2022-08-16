from collections import defaultdict
from functools import partial

import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
import pytransit.qc_display as qc_display
from pytransit.core_data import SessionData, universal
import pytransit

method_wrap_width = 250
selected_export_menu_item = None
convert_menu_item = None
documentation_url = "http://saclab.tamu.edu/essentiality/transit/transit.html"

# sets:
    # universal.selected_method

def create_menu(frame):
    # must imported inside the function to avoid circular import
    import pytransit.components.parameter_panel as parameter_panel
    from pytransit.analysis   import methods as analysis_methods
    from pytransit.export     import methods as export_methods
    from pytransit.convert    import methods as convert_methods
    from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub
    
    global selected_export_menu_item
    global convert_menu_item
    menu_bar = wx.MenuBar(0)
    
    # 
    # File Menu
    # 
    if True:
        file_menu = wx.Menu()
    
        # 
        # Export Menu
        # 
        if True:
            export_menu_item = wx.Menu()
            
            # 
            # Selected Samples
            # 
            if True:
                selected_export_menu_item = wx.Menu()
                export_menu_item.AppendSubMenu(
                    selected_export_menu_item, "Selected Samples"
                )
                
                # 
                # find export options
                # 
                def when_export_clicked(selected_name, event=None):
                    with gui_tools.nice_error_log:
                        if frame.verbose: transit_tools.log(f"Selected Export Method: {selected_name}")
                        gui_tools.run_method_by_label(method_options=export_methods, method_label=selected_name)
                
                for name in export_methods:
                    method = export_methods[name]
                    method.gui.defineMenuItem(frame, method.label)
                    temp_menu_item = method.gui.menuitem
                    selected_export_menu_item.Append(temp_menu_item)
                    
                    frame.Bind(
                        wx.EVT_MENU,
                        partial(when_export_clicked, method.label),
                        temp_menu_item,
                    )

            file_menu.AppendSubMenu(export_menu_item, "Export")

        # 
        # Convert
        # 
        if True:
            convert_menu_item = wx.Menu()
            
            # 
            # prot_table to PTT
            # 
            annotation_convert_pt_to_ptt_menu = wx.MenuItem(
                convert_menu_item,
                wx.ID_ANY,
                "prot_table to PTT",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            convert_menu_item.Append(annotation_convert_pt_to_ptt_menu)
            frame.Bind(wx.EVT_MENU, frame.annotationPT_to_PTT , id=annotation_convert_pt_to_ptt_menu.GetId(),  )
            
            # 
            # prot_table to GFF3
            # 
            annotation_convert_pt_to_gff3_menu = wx.MenuItem(
                convert_menu_item,
                wx.ID_ANY,
                "prot_table to GFF3",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            convert_menu_item.Append(annotation_convert_pt_to_gff3_menu)
            frame.Bind(wx.EVT_MENU, frame.annotationPT_to_GFF3, id=annotation_convert_pt_to_gff3_menu.GetId(), )

            # 
            # PTT to prot_table
            # 
            annotation_convert_ptt_to_pt = wx.MenuItem(
                convert_menu_item,
                wx.ID_ANY,
                "PTT to prot_table",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            convert_menu_item.Append(annotation_convert_ptt_to_pt)
            frame.Bind(wx.EVT_MENU, frame.annotationPTT_to_PT , id=annotation_convert_ptt_to_pt.GetId(),       )
            
            
            # 
            # find Convert options
            # 
            def when_convert_clicked(selected_name, event=None):
                with gui_tools.nice_error_log:
                    if frame.verbose: transit_tools.log(f"Selected Convert Method: {selected_name}")
                    gui_tools.run_method_by_label(method_options=convert_methods, method_label=selected_name)

            for name in convert_methods:
                convert_methods[name].gui.defineMenuItem(frame, convert_methods[name].label)
                temp_menu_item = convert_methods[name].gui.menuitem
                convert_menu_item.Append(temp_menu_item)

                frame.Bind(
                    wx.EVT_MENU,
                    partial(when_convert_clicked, convert_methods[name].label),
                    temp_menu_item,
                )
            
            
            file_menu.AppendSubMenu(convert_menu_item, "Convert")
        
        # 
        # Exit
        # 
        if True:
            exit_option = wx.MenuItem( file_menu, wx.ID_ANY, "&Exit", wx.EmptyString, wx.ITEM_NORMAL )
            file_menu.Append(exit_option)
            def exit_frame(event):
                if frame.verbose: transit_tools.log("Exiting Transit")
                frame.Close()
            frame.Bind(wx.EVT_MENU, exit_frame, id=exit_option.GetId())
        
        menu_bar.Append(file_menu, "&File")
    
    # 
    # View Menu
    # 
    if True:
        view_menu_item = wx.Menu()
        
        # 
        # Scatter Plot
        # 
        if True:
            scatter_menu_item = wx.MenuItem(
                view_menu_item,
                wx.ID_ANY,
                "&Scatter Plot",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            view_menu_item.Append(scatter_menu_item)
            def when_scatter_plot_clicked(event):
                with gui_tools.nice_error_log:
                    from pytransit.components.samples_area import sample_table
                    selected_rows = sample_table.selected_rows
                    if len(dataset_paths) == 2:
                        if frame.verbose: transit_tools.log( f"Showing scatter plot for: {[ each_row['name'] for each_row in selected_rows ]}")
                        data, position = tnseq_tools.CombinedWig.gather_wig_data(
                            list_of_paths=[ each_row["path"] for each_row in selected_rows ],
                        )
                        x = data[0, :]
                        y = data[1, :]

                        plt.plot(x, y, "bo")
                        plt.title("Scatter plot - Reads at TA sites")
                        plt.xlabel(selected_rows[0]["name"])
                        plt.ylabel(selected_rows[1]["name"])
                        plt.show()
                    else:
                        transit_tools.show_error_dialog("Please make sure only two datasets are selected (across control and experimental datasets).")
            frame.Bind(wx.EVT_MENU, when_scatter_plot_clicked, id=scatter_menu_item.GetId() )

        # 
        # Track View
        # 
        if True:
            track_menu_item = wx.MenuItem(view_menu_item, wx.ID_ANY, "&Track View", wx.EmptyString, wx.ITEM_NORMAL)
            view_menu_item.Append(track_menu_item)
            frame.Bind(wx.EVT_MENU, frame.allViewFunc, id=track_menu_item.GetId())
        
        # 
        # Quality Control
        # 
        if True:
            quality_control_option = wx.MenuItem(view_menu_item, wx.ID_ANY, "&Quality Control", wx.EmptyString, wx.ITEM_NORMAL )
            view_menu_item.Append( quality_control_option )
            def when_quality_control_clicked(event):
                with gui_tools.nice_error_log:
                    from pytransit.components.samples_area import sample_table
                    datasets = sample_table.selected_rows # list of dictionaries
                    number_of_files = len(datasets)

                    if number_of_files <= 0:
                        raise Exception(f'''No Datasets selected, unable to run''')
                    else:
                        transit_tools.log(f"Displaying results: {datasets}")
                        try:
                            qc_window = qc_display.qcFrame(frame, datasets)
                            qc_window.Show()
                        except Exception as error:
                            raise Exception(f"Error occured displaying file: {error}")
                        
            frame.Bind(wx.EVT_MENU, when_quality_control_clicked, id=quality_control_option.GetId())
        
        menu_bar.Append(view_menu_item, "&View")
    
    # 
    # Analysis Menu
    # 
    if True:
        analysis_menu = wx.Menu()
        
        # 
        # Himar1 & Tn5
        # 
        if True:
            himar1_menu = wx.Menu()
            tn5_menu = wx.Menu()

            # 
            # generate methods
            # 
            method_names = sorted(analysis_methods.keys())
            for name in method_names:
                method = analysis_methods[name]
                
                def create_callback():
                    # these vars need to be defined here because of how python scopes variables
                    the_method = analysis_methods[name]
                    the_full_name = analysis_methods[name].full_name
                    def load_method_wrapper(event):
                        universal.selected_method = the_method
                        # hide all the other panel stuff
                        for each_method_name in method_names:
                            each_method = analysis_methods[each_method_name]
                            if each_method.gui.panel:
                                each_method.gui.panel.Hide()
                        try:
                            the_method.gui.define_panel(frame)
                        except Exception as error:
                            transit_tools.log(f"Tried to define panel for {the_full_name} but it failed: {error}")
                        return method_select_func(the_full_name, event)
                    return load_method_wrapper
                
                menu_callback = create_callback()
                
                # 
                # himar1 and tn5 menu children
                # 
                for method_name, parent_menu in [ ["himar1", himar1_menu], ["tn5", tn5_menu] ]:
                    temp_menu_item = wx.MenuItem(parent_menu, wx.ID_ANY, method.full_name, wx.EmptyString, wx.ITEM_NORMAL)
                    frame.Bind(wx.EVT_MENU, menu_callback, temp_menu_item)
                    parent_menu.Append(temp_menu_item)
            
            analysis_menu.AppendSubMenu(himar1_menu, "&Himar1 Methods")
            analysis_menu.AppendSubMenu(tn5_menu, "&Tn5 Methods")
        
        menu_bar.Append(analysis_menu, "&Analysis")

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
                    import pytransit.images as images
                    
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


def method_select_func(selected_name, event):
    import pytransit.components.parameter_panel as parameter_panel
    from pytransit.components.parameter_panel import panel
    
    global method_wrap_width
    
    frame = universal.frame
    parameter_panel.hide_all_options()
    
    # If empty is selected
    if selected_name == "[Choose Method]":
        panel.method_info_text.SetLabel("Instructions")
        panel.method_instructions.Show()
        panel.method_instructions.SetLabel(frame.instructions_text)
        panel.method_instructions.Wrap(method_wrap_width)
        panel.method_short_text.Hide()
        panel.method_long_text.Hide()
        panel.method_tn_text.Hide()
        panel.method_desc_text.Hide()

        panel.method_choice = ""
    else:
        panel.method_sizer_text.Show()
        
        from pytransit.analysis import methods as analysis_methods
        
        matched_name = None
        # Get selected Method and hide Others
        for name in analysis_methods:
            try: analysis_methods[name].gui.panel.Hide()
            except Exception as error: pass
            
            if analysis_methods[name].full_name == selected_name:
                matched_name = name
        
        if matched_name in analysis_methods:
            name = matched_name
            panel.method_info_text.SetLabel("%s" % analysis_methods[name].long_name)

            panel.method_tn_text.Show()
            panel.method_tn_text.SetLabel(analysis_methods[name].transposons_text)
            panel.method_tn_text.Wrap(method_wrap_width)

            panel.method_desc_text.Show()
            panel.method_desc_text.SetLabel(analysis_methods[name].long_desc)
            panel.method_desc_text.Wrap(method_wrap_width)
            panel.method_instructions.SetLabel(" ")
            analysis_methods[name].gui.panel.Show()
            frame.statusBar.SetStatusText("[%s]" % analysis_methods[name].short_name)

        parameter_panel.show_progress_section()
        panel.method_choice = selected_name

    frame.Layout()
    if frame.verbose:
        transit_tools.log("Selected Method: %s" % (selected_name))
