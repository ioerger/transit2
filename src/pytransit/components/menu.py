from collections import defaultdict
from functools import partial

import pytransit.gui_tools as gui_tools
import pytransit.transit_tools as transit_tools
import pytransit.qc_display as qc_display
from pytransit.core_data import SessionData, universal

selected_export_menu_item = None
convert_menu_item = None

# sets:
    # universal.selected_method

def create_menu(frame):
    # must be done here to avoid circular import
    import pytransit.components.parameter_panel as parameter_panel
    from pytransit.analysis   import methods
    from pytransit.export     import methods as export_methods
    from pytransit.convert    import methods as convert_methods
    from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub
    
    global selected_export_menu_item
    global convert_menu_item
    
    # 
    # define visual elements
    # 
    if True:
        menu_bar         = wx.MenuBar(0)
        
        # 
        # File Menu
        # 
        if True:
            file_menu   = wx.Menu()
        
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
                
                
                # 
                # find Convert options
                # 
                def when_convert_clicked(self, selected_name, event=None):
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
                
                
                convert_menu_item.Append(annotation_convert_ptt_to_pt)
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
            scatter_menu_item = wx.MenuItem(
                view_menu_item,
                wx.ID_ANY,
                "&Scatter Plot",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )

            view_menu_item.Append(scatter_menu_item)

            track_menu_item = wx.MenuItem(
                view_menu_item, wx.ID_ANY, "&Track View", wx.EmptyString, wx.ITEM_NORMAL
            )

            view_menu_item.Append(track_menu_item)
            menu_bar.Append(view_menu_item, "&View")
        
        # 
        # Methods Menu
        # 
        if True:
            methods_menu_item = wx.Menu()
            himar1_menu_item = wx.Menu()
            tn5_menu_item = wx.Menu()

            methods_menu_item.AppendSubMenu(himar1_menu_item, "&Himar1 Methods")
            methods_menu_item.AppendSubMenu(tn5_menu_item, "&Tn5 Methods")
            menu_bar.Append(methods_menu_item, "&Analysis")

        # 
        # Help Menu
        # 
        if True:
            help_menu_item = wx.Menu()
            documentation_menu_item = wx.MenuItem(
                help_menu_item,
                wx.ID_ANY,
                "&Documentation",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            help_menu_item.Append(documentation_menu_item)
            about_menu_item = wx.MenuItem(
                help_menu_item, wx.ID_ANY, "&About", wx.EmptyString, wx.ITEM_NORMAL
            )
            help_menu_item.Append(about_menu_item)

            menu_bar.Append(help_menu_item, "&Help")
        
        frame.SetMenuBar(menu_bar)
    
    # 
    # events
    # 
    if True:
        frame.Bind(wx.EVT_MENU, frame.annotationPT_to_PTT , id=annotation_convert_pt_to_ptt_menu.GetId(),  )
        frame.Bind(wx.EVT_MENU, frame.annotationPT_to_GFF3, id=annotation_convert_pt_to_gff3_menu.GetId(), )
        frame.Bind(wx.EVT_MENU, frame.annotationPTT_to_PT , id=annotation_convert_ptt_to_pt.GetId(),       )
        frame.Bind(wx.EVT_MENU, frame.scatterFunc         , id=scatter_menu_item.GetId()                   )
        frame.Bind(wx.EVT_MENU, frame.allViewFunc         , id=track_menu_item.GetId()                     )
        frame.Bind(wx.EVT_MENU, frame.aboutFunc           , id=about_menu_item.GetId()                  )
        frame.Bind(wx.EVT_MENU, frame.documentationFunc   , id=documentation_menu_item.GetId()          )
    
    # 
    # generate methods
    # 
    method_names = sorted(methods.keys())
    for name in method_names:
        method = methods[name]
        
        def create_callback():
            # these vars need to be defined here because of how python scopes variables
            the_method = methods[name]
            the_full_name = methods[name].full_name
            def load_method_wrapper(event):
                universal.selected_method = the_method
                # hide all the other panel stuff
                for each_method_name in method_names:
                    each_method = methods[each_method_name]
                    if each_method.gui.panel:
                        each_method.gui.panel.Hide()
                try:
                    the_method.gui.define_panel(frame)
                except Exception as error:
                    pass
                return method_select_func(the_full_name, event)
            return load_method_wrapper
        
        menu_callback = create_callback()
        
        # attach menus
        for method_name, parent_menu in [ ["himar1", himar1_menu_item], ["tn5", tn5_menu_item] ]:
            temp_menu_item = wx.MenuItem(parent_menu, wx.ID_ANY, method.full_name, wx.EmptyString, wx.ITEM_NORMAL,)
            frame.Bind(wx.EVT_MENU,menu_callback,temp_menu_item,)
            parent_menu.Append(temp_menu_item)


def method_select_func(selected_name, event):
    import pytransit.components.parameter_panel as parameter_panel
    from pytransit.components.parameter_panel import panel
    
    method_wrap_width = 250
    frame = universal.frame
    parameter_panel.hide_all_options()
    
    # If empty is selected
    if selected_name == "[Choose Method]":
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
        panel.method_sizer_text.Show()
        
        from pytransit.analysis import methods
        
        matched_name = None
        # Get selected Method and hide Others
        for name in methods:
            try: methods[name].gui.panel.Hide()
            except Exception as error: pass
            
            if methods[name].full_name == selected_name:
                matched_name = name
        
        if matched_name in methods:
            name = matched_name
            panel.method_info_text.SetLabel("%s" % methods[name].long_name)

            panel.method_tn_text.Show()
            panel.method_tn_text.SetLabel(methods[name].transposons_text)
            panel.method_tn_text.Wrap(method_wrap_width)

            panel.method_desc_text.Show()
            panel.method_desc_text.SetLabel(methods[name].long_desc)
            panel.method_desc_text.Wrap(method_wrap_width)
            panel.method_instructions.SetLabel(" ")
            methods[name].gui.panel.Show()
            frame.statusBar.SetStatusText("[%s]" % methods[name].short_name)

        parameter_panel.show_progress_section()
        panel.method_choice = selected_name

    frame.Layout()
    if frame.verbose:
        transit_tools.log("Selected Method: %s" % (selected_name))


def method_select_helper(method, event):
    import pytransit.components.parameter_panel as parameter_panel
    from pytransit.components.parameter_panel import panel
    
    frame = universal.frame
    
    panel.method_sizer_text.Show()
    
    from pytransit.analysis import methods
    
    # Hide all the others
    for each_key, each in var.items():
        if each.panel:
            each.panel.Hide()
    
    panel.method_info_text.SetLabel(method.long_name)

    panel.method_tn_text.Show()
    panel.method_tn_text.SetLabel(method.transposons_text)
    panel.method_tn_text.Wrap(250)

    panel.method_desc_text.Show()
    panel.method_desc_text.SetLabel(method.long_desc)
    panel.method_desc_text.Wrap(250)
    panel.method_instructions.SetLabel(" ")
    
    method.define_panel()
    
    frame.statusBar.SetStatusText(f"[{method.short_name}]")

    parameter_panel.show_progress_section()
    panel.method_choice = selected_name

    frame.Layout()
    if frame.verbose:
        transit_tools.log("Selected Method: %s" % (selected_name))

