from collections import defaultdict
from functools import partial

selected_export_menu_item = None
convert_menu_item = None

def create_menu(self):
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
        file_menu_item   = wx.Menu()
        export_menu_item = wx.Menu()
        selected_export_menu_item = wx.Menu()

        # Selected datasets
        export_menu_item.AppendSubMenu(
            selected_export_menu_item, "Selected Datasets"
        )

        file_menu_item.AppendSubMenu(export_menu_item, "Export")

        convert_menu_item = wx.Menu()
        annotation_convert_pt_to_ptt_menu = wx.MenuItem(
            convert_menu_item,
            wx.ID_ANY,
            "prot_table to PTT",
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )
        convert_menu_item.Append(annotation_convert_pt_to_ptt_menu)

        annotation_convert_pt_to_gff3_menu = wx.MenuItem(
            convert_menu_item,
            wx.ID_ANY,
            "prot_table to GFF3",
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )
        convert_menu_item.Append(annotation_convert_pt_to_gff3_menu)

        annotation_convert_ptt_to_pt = wx.MenuItem(
            convert_menu_item,
            wx.ID_ANY,
            "PTT to prot_table",
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )

        convert_menu_item.Append(annotation_convert_ptt_to_pt)

        # self.annotationConvertGFF3ToPT = wx.MenuItem( convert_menu_item, wx.ID_ANY, "GFF3 to prot_table", wx.EmptyString, wx.ITEM_NORMAL )
        # convert_menu_item.Append( self.annotationConvertGFF3ToPT )
        file_menu_item.AppendSubMenu(convert_menu_item, "Convert")

        file_exit_menu_item = wx.MenuItem(
            file_menu_item, wx.ID_ANY, "&Exit", wx.EmptyString, wx.ITEM_NORMAL
        )
        file_menu_item.Append(file_exit_menu_item)
        menu_bar.Append(file_menu_item, "&File")

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
        methods_menu_item = wx.Menu()
        himar1_menu_item = wx.Menu()
        tn5_menu_item = wx.Menu()

        methods_menu_item.AppendSubMenu(himar1_menu_item, "&Himar1 Methods")
        methods_menu_item.AppendSubMenu(tn5_menu_item, "&Tn5 Methods")
        menu_bar.Append(methods_menu_item, "&Analysis")

        self.SetMenuBar(menu_bar)

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

    # 
    # Export options
    # 
    for name in export_methods:
        export_methods[name].gui.defineMenuItem(self, export_methods[name].label)
        temp_menu_item = export_methods[name].gui.menuitem
        selected_export_menu_item.Append(temp_menu_item)

        self.Bind(
            wx.EVT_MENU,
            partial(self.ExportSelectFunc, export_methods[name].label),
            temp_menu_item,
        )
    
    # 
    # Convert options
    # 
    for name in convert_methods:
        convert_methods[name].gui.defineMenuItem(self, convert_methods[name].label)
        temp_menu_item = convert_methods[name].gui.menuitem
        convert_menu_item.Append(temp_menu_item)

        self.Bind(
            wx.EVT_MENU,
            partial(self.ConvertSelectFunc, convert_methods[name].label),
            temp_menu_item,
        )
    
    # 
    # events
    # 
    if True:
        self.Bind(wx.EVT_MENU, self.annotationPT_to_PTT , id=annotation_convert_pt_to_ptt_menu.GetId(),  )
        self.Bind(wx.EVT_MENU, self.annotationPT_to_GFF3, id=annotation_convert_pt_to_gff3_menu.GetId(), )
        self.Bind(wx.EVT_MENU, self.annotationPTT_to_PT , id=annotation_convert_ptt_to_pt.GetId(),       )
        self.Bind(wx.EVT_MENU, self.Exit                , id=file_exit_menu_item.GetId()                 )
        self.Bind(wx.EVT_MENU, self.scatterFunc         , id=scatter_menu_item.GetId()                   )
        self.Bind(wx.EVT_MENU, self.allViewFunc         , id=track_menu_item.GetId()                     )
        self.Bind(wx.EVT_MENU, self.aboutFunc           , id=about_menu_item.GetId()                  )
        self.Bind(wx.EVT_MENU, self.documentationFunc   , id=documentation_menu_item.GetId()          )
    
    
    method_order = [("gumbel", 1), ("resampling", 2), ("hmm", 3)]
    order = defaultdict(lambda: 100)
    for k, v in method_order:
        order[k] = v

    for name in sorted(methods.keys(), key=lambda x: order[x]):
        methods[name].gui.define_panel(self)
        methods[name].gui.Hide()

        if "himar1" in methods[name].transposons:
            temp_menu_item = wx.MenuItem(
                himar1_menu_item,
                wx.ID_ANY,
                methods[name].fullname(),
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            self.Bind(
                wx.EVT_MENU,
                partial(method_select_func, methods[name].fullname()),
                temp_menu_item,
            )

            himar1_menu_item.Append(temp_menu_item)

        if "tn5" in methods[name].transposons:
            temp_menu_item = wx.MenuItem(
                tn5_menu_item,
                wx.ID_ANY,
                methods[name].fullname(),
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            self.Bind(
                wx.EVT_MENU,
                partial(method_select_func, methods[name].fullname()),
                temp_menu_item,
            )
            tn5_menu_item.Append(temp_menu_item)

def method_select_func(selected_name, event):
    import pytransit.components.parameter_panel as parameter_panel
    from pytransit.components.parameter_panel import panel
    
    # If empty is selected
    if selected_name == "[Choose Method]":
        method_wrap_width = 250
        parameter_panel.hide_all_options()
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
        parameter_panel.show_global_options()
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

        parameter_panel.show_progress_section()
        panel.method_choice = selected_name

    frame.Layout()
    if frame.verbose:
        transit_tools.transit_message("Selected Method: %s" % (selected_name))

