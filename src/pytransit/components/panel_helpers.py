from pytransit.specific_tools.transit_tools import wx, GenBitmapTextButton
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools import gui_tools, transit_tools, tnseq_tools, norm_tools, stat_tools

default_padding = 30
default_label_size = (220, -1)
default_widget_size = (100, -1)

# 
# 
# generic
# 
# 
if True:
    class NewPanel:
        def __init__(self, *args, **kwargs):
            pass
        
        def __enter__(self):
            self.wx_panel = wx.lib.scrolledpanel.ScrolledPanel(
                gui.frame,
                wx.ID_ANY,
                wx.DefaultPosition,
                # wx.DefaultSize,
                wx.Size(int(gui.width/4), int(gui.height*0.35)),
                wx.TAB_TRAVERSAL,
            )
            # self.wx_panel.SetMaxSize((width + (width - width_2) + dx, -1)) # Trying to limit the width of our frame
            self.main_sizer = wx.BoxSizer(wx.VERTICAL)
            NewPanel.recent = self
            return (self.wx_panel, self.main_sizer)
        
        def __exit__(self, _, error, traceback_obj):
            if error is not None:
                import traceback
                print(''.join(traceback.format_tb(traceback_obj)))
                frame = gui.frame
                if frame and hasattr(frame, "status_bar"):
                    frame.status_bar.SetStatusText("Error: "+str(error.args))
            else:
                self.refresh()
        
        def refresh(self):
            from pytransit.components import parameter_panel
            parameter_panel.set_panel(self.wx_panel)
            self.wx_panel.SetSizer(self.main_sizer)
            self.wx_panel.Layout()
            self.main_sizer.Fit(self.wx_panel)
            self.wx_panel.SetupScrolling()
            gui.frame.Layout()
            
    
    def create_tooltip_and_label(panel, tooltip_text, label_text=None):
        from pytransit.components.icon import InfoIcon
        inner_sizer = wx.BoxSizer(wx.HORIZONTAL)
        inner_sizer.Add(
            InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
            proportion=0,
            flag=wx.ALIGN_CENTER_VERTICAL,
            border=gui_tools.default_padding,
        )
        # a spacer because the border doesn't seem to actually work
        inner_sizer.Add(10, default_padding)
        if label_text != None:
            label = wx.StaticText(panel, wx.ID_ANY, label_text, wx.DefaultPosition, default_label_size, 0)
            label.Wrap(-1)
            inner_sizer.Add(
                label,
                proportion=0,
                flag=wx.ALIGN_CENTER_VERTICAL,
                border=gui_tools.default_padding,
            )
        return inner_sizer
    
    def create_button(panel, sizer, *, label):
        """
        Example:
            @panel_helpers.create_button(pop_up_panel, main_sizer, label="")
            def when_button_clicked(event):
                print("do stuff")
        """
        a_button = wx.Button(
            panel,
            wx.ID_ANY,
            label,
            wx.DefaultPosition,
            wx.DefaultSize,
            0,
        )
        sizer.Add(a_button, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL, border=gui_tools.default_padding)
        def decorator(func):
            @gui_tools.bind_to(a_button, wx.EVT_BUTTON)
            def wrapper(*args,**kwargs):
                with gui_tools.nice_error_log:
                    return func(*args,**kwargs)
            return wrapper
        return decorator

    

    def create_generic_button_to_popup(panel, sizer, *, button_label, tooltip_text="", popup_title=""):
        """
        On click of a button, show popup with a dropdown and ok button
        On "OK", it will return the value you pass in
        """

        from os.path import basename, dirname
        row_sizer = create_tooltip_and_label(panel, tooltip_text=tooltip_text)
        if True:
            # the button to click on in passed in panel to get a popup
            click_button = wx.Button(
                panel,
                wx.ID_ANY,
                button_label,
                wx.DefaultPosition,
                wx.DefaultSize,
                0,
            )
            results=None
            @gui_tools.bind_to(click_button, wx.EVT_BUTTON)
            def when_button_clicked(*args,**kwargs):
                nonlocal results
                win = wx.Dialog(panel,wx.FRAME_FLOAT_ON_PARENT)
                popup_sizer = wx.BoxSizer(wx.VERTICAL)
                win.SetSizer(popup_sizer)

                dropdown_label_text= wx.StaticText(win, wx.ID_ANY, label="Select One : ", style=wx.ALIGN_LEFT)
                popup_sizer.Add(dropdown_label_text, 0, wx.ALL, gui_tools.default_padding)
                dropdown_selection = wx.ComboBox(win,choices = ["Yes", "No", "Maybe"])
                popup_sizer.Add(dropdown_selection, wx.ALL, gui_tools.default_padding)

                ok_btn = wx.Button(win, wx.ID_OK, label = "Ok", size = (50,20), pos = (75,50))
                popup_sizer.Add(ok_btn,wx.EXPAND, gui_tools.default_padding)

                win.Layout()
                popup_sizer.Fit(win)
                res = win.ShowModal()
                if res == wx.ID_OK:
                    results = dropdown_selection.GetValue()
                win.Destroy()
        row_sizer.Add(click_button, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
        
        sizer.Add(row_sizer, 0, wx.ALIGN_LEFT, gui_tools.default_padding)
        return lambda *args, **kwargs: results


    def create_folder_input(panel, sizer, *, button_label, init_folder_text="", tooltip_text="", popup_title="", default_folder=None, default_folder_name="", allowed_extensions='All folders (*.*)|*.*', after_select=lambda *args: None, alignment=wx.ALIGN_CENTER, size=(-1,-1), color=gui_tools.color.light_gray):
        """
            Example:
                folder_path_getter = create_folder_input(self.panel, main_sizer, button_label="Add context folder", allowed_extensions='All folders (*.*)|*.*')
                folder_path_or_none = folder_path_getter()
        """
        from os.path import basename
        row_sizer = create_tooltip_and_label(panel, tooltip_text=tooltip_text)
        if True:
            # 
            # button
            # 
            if True:
                add_folder_button = GenBitmapTextButton(
                    panel,
                    wx.ID_ANY,
                    gui_tools.bit_map,
                    button_label,
                    size=wx.Size(*size),
                )
                add_folder_button.SetMinSize(size)
                add_folder_button.SetMaxSize(size)
                add_folder_button.SetBackgroundColour(color)
                
                folder_text = None
                the_folder_path = None
                # whenever the button is clicked, set the folder
                @gui_tools.bind_to(add_folder_button, wx.EVT_BUTTON)
                def when_button_clicked(*args,**kwargs):
                    nonlocal the_folder_path
                    with gui_tools.nice_error_log:
                        # set the folder path variable
                        the_folder_path = gui_tools.ask_for_folder(
                            message=popup_title,
                            default_folder=default_folder,
                        )
                        folder_text.SetLabel(basename(the_folder_path or ""))
                        after_select(*args, the_folder_path)
            row_sizer.Add(add_folder_button, 0, wx.ALIGN_CENTER, gui_tools.default_padding)
            # padding
            row_sizer.Add(20, 1)
            
            # 
            # Text
            # 
            folder_text = wx.StaticText(panel, wx.ID_ANY, label=init_folder_text, style=wx.ALIGN_LEFT)
            row_sizer.Add(folder_text, 0, wx.ALIGN_CENTER, gui_tools.default_padding)
        
        sizer.Add(row_sizer, proportion=0, flag=alignment, border=gui_tools.default_padding)
        output = lambda *args, **kwargs: the_folder_path
        output.set_label = lambda text: folder_text.SetLabel(text)
        return output
    
    def create_file_input(panel, sizer, *, button_label, init_file_text="", tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*', after_select=lambda *args: None, alignment=wx.ALIGN_CENTER, size=(-1,-1), color=gui_tools.color.light_gray):
        """
            Example:
                file_path_getter = create_file_input(self.panel, main_sizer, button_label="Add context file", allowed_extensions='All files (*.*)|*.*')
                file_path_or_none = file_path_getter()
        """
        from os.path import basename, dirname
        row_sizer = create_tooltip_and_label(panel, tooltip_text=tooltip_text)
        if True:
            # 
            # button
            # 
            if True:
                add_file_button = GenBitmapTextButton(
                    panel,
                    wx.ID_ANY,
                    gui_tools.bit_map,
                    button_label,
                    size=wx.Size(*size),
                )
                add_file_button.SetMinSize(size)
                add_file_button.SetMaxSize(size)
                add_file_button.SetBackgroundColour(color)
                
                file_text = None
                the_file_path = None
                # whenever the button is clicked, set the file
                @gui_tools.bind_to(add_file_button, wx.EVT_BUTTON)
                def when_button_clicked(*args,**kwargs):
                    nonlocal the_file_path
                    with gui_tools.nice_error_log:
                        # set the file path variable
                        the_file_path = gui_tools.ask_for_file(
                            message=popup_title,
                            default_folder=default_folder,
                            default_file_name=default_file_name,
                            allowed_extensions=allowed_extensions,
                        )
                        name = basename(the_file_path or "")
                        parent_name = dirname(the_file_path or "")
                        if parent_name != "." and parent_name != "":
                            name = basename(parent_name)+"/"+name
                        file_text.SetLabel(name)
                        after_select(*args, the_file_path)
            row_sizer.Add(add_file_button, 0, wx.ALIGN_CENTER, gui_tools.default_padding)
            # padding
            row_sizer.Add(20, 1)
            
            # 
            # Text
            # 
            file_text = wx.StaticText(panel, wx.ID_ANY, label=init_file_text, style=wx.ALIGN_LEFT)
            row_sizer.Add(file_text, 0, wx.ALIGN_CENTER, gui_tools.default_padding)
        
        sizer.Add(row_sizer, proportion=0, flag=alignment, border=gui_tools.default_padding)
        output = lambda *args, **kwargs: the_file_path
        output.set_label = lambda text: file_text.SetLabel(text)
        return output
    
    def create_multi_file_input(panel, sizer, *, button_label, tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*', after_select=lambda *args: None):
        """
            Example:
                file_path_getter = create_file_input(self.panel, main_sizer, button_label="Add context file", allowed_extensions='All files (*.*)|*.*')
                file_path_or_none = file_path_getter()
        """
        from os.path import basename
        row_sizer = create_tooltip_and_label(panel, tooltip_text=tooltip_text)
        if True:
            # 
            # button
            # 
            if True:
                add_file_button = wx.Button(
                    panel,
                    wx.ID_ANY,
                    button_label,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                file_text = None
                selected_paths = None
                # whenever the button is clicked, set the file
                @gui_tools.bind_to(add_file_button, wx.EVT_BUTTON)
                def when_button_clicked(*args,**kwargs):
                    nonlocal selected_paths
                    with gui_tools.nice_error_log:
                        # set the file path variable
                        selected_paths = gui_tools.ask_for_files(
                            message=popup_title,
                            default_folder=default_folder,
                            default_file_name=default_file_name,
                            allowed_extensions=allowed_extensions,
                        )
                        names = "\n".join([ basename(each) for each in selected_paths ])
                        file_text.SetLabel(names)
                        after_select(*args)
            row_sizer.Add(add_file_button, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
            
            # 
            # Text
            # 
            file_text = wx.StaticText(panel, wx.ID_ANY, label="", style=wx.ALIGN_LEFT)
            row_sizer.Add(file_text, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
        
        sizer.Add(row_sizer, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL, border=gui_tools.default_padding)
        return lambda *args, **kwargs: selected_paths
    
    def create_persistent_file_input(panel, sizer, *, name, button_label, tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*'):
        """
            Example:
                file_path_getter = create_persistent,file_input(self.panel, main_sizer, button_label="Add context file", allowed_extensions='All files (*.*)|*.*')
                file_path_or_none = file_path_getter()
        """
        ##todo : HANDLE NAME
        from os.path import basename
        row_sizer = create_tooltip_and_label(panel, tooltip_text=tooltip_text)
        if True:
            # 
            # button
            # 
            if True:
                add_file_button = wx.Button(
                    panel,
                    wx.ID_ANY,
                    button_label,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                file_text = None
                the_file_path = None
                # whenever the button is clicked, set the file
                @gui_tools.bind_to(add_file_button, wx.EVT_BUTTON)
                def when_button_clicked(*args,**kwargs):
                    nonlocal the_file_path
                    with gui_tools.nice_error_log:
                        # set the file path variable
                        the_file_path = gui_tools.ask_for_file(
                            message=popup_title,
                            default_folder=default_folder,
                            default_file_name=default_file_name,
                            allowed_extensions=allowed_extensions,
                        )
                        file_text.SetLabel(basename(the_file_path or ""))
            row_sizer.Add(add_file_button, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
            
            # 
            # Text
            # 
            file_text = wx.StaticText(panel, wx.ID_ANY, label="", style=wx.ALIGN_LEFT)
            row_sizer.Add(file_text, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
            
        sizer.Add(row_sizer, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL, border=gui_tools.default_padding)
        return lambda *args, **kwargs: the_file_path
   
    def define_choice_box(
        panel,
        *,
        label_text="",
        options=[""],
        tooltip_text="",
        label_size=None,
        widget_size=None,
    ):
        from pytransit.components.icon import InfoIcon
        
        if not label_size:
            label_size = default_label_size
        if not widget_size:
            widget_size = default_widget_size

        sizer = create_tooltip_and_label(panel, tooltip_text=tooltip_text, label_text=label_text)
        choice_box = wx.Choice(panel, wx.ID_ANY, wx.DefaultPosition, widget_size, options, 0 )
        choice_box.SetSelection(0)
        sizer.Add(choice_box, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
        return (None, choice_box, sizer)
    
    def create_label(
        panel,
        sizer,
        text="",
        size=None,
    ):
        if not label_size:
            label_size = default_label_size
        inner_sizer = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(panel, wx.ID_ANY, label_text, wx.DefaultPosition, label_size, 0)
        label.Wrap(-1)
        choice_box = wx.Choice(panel, wx.ID_ANY, wx.DefaultPosition, widget_size, options, 0 )
        choice_box.SetSelection(0)
        inner_sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
        
        sizer.Add(inner_sizer, 1, wx.ALIGN_LEFT, gui_tools.default_padding)

    def create_choice_input(panel, sizer, label, options, default_option=None, tooltip_text=""):
        # 
        # setup arguments
        # 
        if len(options) == 0:
            options.append("")
        if default_option is None:
            default_option = options[0]
        # 
        # create box
        # 
        (
            label,
            wxobj,
            inner_sizer,
        ) = define_choice_box(
            panel,
            label_text=label,
            options=options,
            tooltip_text=tooltip_text,
        )
        sizer.Add(inner_sizer, 1, wx.ALIGN_LEFT, gui_tools.default_padding)
        # return a value-getter
        wxobj.SetSelection(wxobj.FindString(default_option))
        return lambda *args: wxobj.GetString(wxobj.GetCurrentSelection())
    
    def define_text_box(
            panel,
            label_text="",
            default_value="",
            tooltip_text="",
            label_size=None,
            widget_size=None,
        ):
            from pytransit.components.icon import InfoIcon
            if not label_size:
                label_size = default_label_size
            if not widget_size:
                widget_size = default_widget_size

            sizer = create_tooltip_and_label(panel, tooltip_text=tooltip_text, label_text=label_text)
            text_box = wx.TextCtrl(panel, wx.ID_ANY, f"{default_value}", wx.DefaultPosition, widget_size, 0)
            sizer.Add(text_box, 0,  wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
            sizer.Layout()
            
            return None, text_box, sizer
            
    def create_check_box_getter(panel, sizer, *, label_text="", default_value=False, tooltip_text="", label_size=None, widget_size=None):
        from pytransit.components.icon import InfoIcon
        if not widget_size:
            widget_size = default_widget_size
        if not label_size:
            label_size = default_label_size
        
        inner_sizer = create_tooltip_and_label(panel, tooltip_text=tooltip_text, label_text=label_text)
        check_box   = wx.CheckBox(panel, label="", size=widget_size)
        check_box.SetValue(default_value)
        inner_sizer.Add(check_box, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
        
        sizer.Add(inner_sizer, 1, wx.ALIGN_LEFT, gui_tools.default_padding)
        return lambda *args: check_box.GetValue()
    
    def create_text_box_getter(panel, sizer, label_text="", default_value="", tooltip_text="", label_size=None, widget_size=None,):
        (
            _,
            wxobj,
            wrapper_sizer,
        ) = define_text_box(
            panel,
            label_text=label_text,
            default_value=str(default_value),
            tooltip_text=tooltip_text,
        )
        sizer.Add(wrapper_sizer, 1, wx.ALIGN_LEFT, gui_tools.default_padding)
        return lambda *args: wxobj.GetValue()

    def create_float_getter(panel, sizer, *, label_text, default_value, tooltip_text=None):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text=label_text,
            default_value=str(float(default_value)),
            tooltip_text=tooltip_text,
        )
        return lambda *args: float(get_text())        
    
    def create_int_getter(panel, sizer, *, label_text, default_value, tooltip_text=None):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text=label_text,
            default_value=str(int(default_value)),
            tooltip_text=tooltip_text,
        )
        return lambda *args: int(get_text())        
    
    def create_multiselect_getter(panel, sizer, options, label_text=None, tooltip_text=None):
        from pytransit.components import parameter_panel
        from pytransit.components.generic.table import Table
        
        sizer.Add(
            create_tooltip_and_label(panel, tooltip_text=tooltip_text, label_text=label_text),
            0,
            wx.ALIGN_LEFT, gui_tools.default_padding
        )
        
        table = None
        row_height_approximate = 25
        with Table(frame=panel, column_width="100%", min_size=(parameter_panel.panel.max_width*0.5, -1)) as table:
            sizer.Add(
                table.wx_object,
                0,
                wx.ALIGN_CENTER_HORIZONTAL,
                gui_tools.default_padding
            )
            for each in options:
                table.add(dict(option=each))
        
        return lambda *args: [ each["option"] for each in table.selected_rows ]

# 
# 
# specific
# 
# 
if True:
    def create_normalization_input(panel, sizer, default="TTR"):
        from pytransit.methods.normalize import Method
        return create_choice_input(
            panel,
            sizer,
            label="Normalization: ",
            options=Method.options,
            tooltip_text="Choice of normalization method. The default choice is TTR (trimmed total reads). See documentation for a description other methods."
        )
    
    def create_wig_choice(panel, sizer, *, label_text, tooltip_text="Pick any wig inside of the combined wig"):
        samples_at_the_time = gui.samples
        wig_ids = [ each.id for each in samples_at_the_time ]
        getter_of_wig_id = create_choice_input(
            panel,
            sizer,
            label=label_text,
            options=wig_ids,
            tooltip_text=tooltip_text,
        )
        return lambda *args: samples_at_the_time[wig_ids.index(getter_of_wig_id(*args))]
    
    
    def create_reference_condition_input(panel, sizer):
        return create_choice_input(
            panel,
            sizer,
            label="Ref Condition:",
            options=[ "[None]" ] + [ each.name for each in gui.conditions ],
            tooltip_text="Which condition to use as a reference for calculating LFCs. If no ref given, this compares against grand mean of all conditions.",
        )
    
    def create_condition_input(panel, sizer, *, label_text="Condition", tooltip_text="choose condition"):
        return create_choice_input(
            panel,
            sizer,
            label=label_text,
            options=[ "[None]" ] + [ each.name for each in gui.conditions],
            tooltip_text=tooltip_text,
        )
        
    def create_include_condition_list_input(panel, sizer):
        (
            _,
            wxobj,
            wrapper_sizer,
        ) = define_text_box(
            panel=panel,
            label_text="Include Conditions\n",
            default_value="",
            tooltip_text="comma seperated list (default=all)",
        )
        sizer.Add(wrapper_sizer, 1, wx.ALIGN_LEFT, gui_tools.default_padding)
        
        def get_value(*args):
            as_list = wxobj.GetValue().split(",")
            without_empty_strings = [ each for each in as_list if len(each) > 0 ]
            if len(without_empty_strings) == 0:
                return [ each.name for each in gui.conditions ]
            else:
                return without_empty_strings
        
        return get_value
    
    def create_exclude_condition_list_input(panel, sizer):
        (
            _,
            wxobj,
            wrapper_sizer,
        ) = define_text_box(
            panel=panel,
            label_text="Exclude Conditions\n",
            default_value="",
            tooltip_text="comma seperated list (default=none)",
        )
        sizer.Add(wrapper_sizer, 1, wx.ALIGN_LEFT, gui_tools.default_padding)
        
        def get_value(*args):
            as_list = wxobj.GetValue().split(",")
            without_empty_strings = [ each for each in as_list if len(each) > 0 ]
            return without_empty_strings
        
        return get_value
    
    def create_control_condition_input(panel, sizer):
        (
            label,
            ref_condition_wxobj,
            ref_condition_choice_sizer,
        ) = define_choice_box(
            panel,
            label_text="Control Condition:",
            options=[ "[None]" ] + [ each.name for each in gui.conditions ],
            tooltip_text="In a comparison of a treatment vs. control condition, this is the control condition to be used as reference",
            label_size=(200, 20),
        )
        sizer.Add(ref_condition_choice_sizer, 1, wx.ALIGN_LEFT, gui_tools.default_padding)
        return lambda *args: ref_condition_wxobj.GetString(ref_condition_wxobj.GetCurrentSelection())
    
    def create_experimental_condition_input(panel, sizer):
        (
            label,
            ref_condition_wxobj,
            ref_condition_choice_sizer,
        ) = define_choice_box(
            panel,
            label_text="Experimental Condition:",
            options=[ "[None]" ] + [ each.name for each in gui.conditions ],
            tooltip_text="In a comparison of a treatment vs. control condition, this is the treatment condition”p",
            label_size=(200, 20),
        )
        sizer.Add(ref_condition_choice_sizer, 1, wx.ALIGN_LEFT, gui_tools.default_padding)
        return lambda *args: ref_condition_wxobj.GetString(ref_condition_wxobj.GetCurrentSelection())
    
    def create_n_terminus_input(panel, sizer):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text="Ignore N-Terminus %:",
            default_value="0",
            tooltip_text="Ignore TA sites in a given fraction of the N-terminal end of the ORF. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",    
        )
        return lambda *args: float(get_text())
    
    def create_c_terminus_input(panel, sizer):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text="Ignore C-Terminus %:",
            default_value="0",
            tooltip_text="Ignore TA sites in a given fraction of the C-terminal end of the ORF. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",    
        )
        return lambda *args: float(get_text())
    
    def create_pseudocount_input(panel, sizer, default_value="5", tooltip = "Pseudo-counts used in calculating log-fold-change. Note: pseudocounts do not affect P values. Useful to dampen the effects of small counts which may lead to deceptively high LFC."):
        # 
        # text input: Pseudocount
        # 
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text="Pseudocount:",
            default_value=default_value,
            tooltip_text=tooltip
        )
        return lambda *args: float(get_text())
    
    
    def create_alpha_input(panel, sizer, default_value="1000"):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text="Alpha:",
            default_value=default_value,
            tooltip_text="Value added to MSE in F-test for moderated ANOVA: F = MSR/(MSE+alpha). Note: this value does affect P values. This is helpful because genes with very low counts are occasionally ranked as significant by traditional ANOVA, even though the apparent variability is probably due to noise. Setting alpha to a number like 1000 helps filter out these irrelevant genes by reducing their significance. If you want to emulate the standard ANOVA test, you can set alpha to 0.",    
        )
        return lambda *args: int(get_text())
    
    def create_winsorize_input(panel, sizer, default_value=False):
        return create_check_box_getter(
            panel,
            sizer,
            label_text="Winsorize:",
            default_value=default_value,
            tooltip_text="Winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers).",    
        )
    
    def create_selected_condition_names_input(panel, sizer, default_value=False):
        check_box_getter = create_check_box_getter(panel, sizer,
            label_text="Use only selected conditions",
            default_value=default_value,
            tooltip_text="When checked, use the conditions table (on the left) to select which conditions to run this analysis on. When not checked, all conditions are used.",
        )
        def wrapper(*args, **kwargs):
            is_checked = check_box_getter(*args, **kwargs)
            if is_checked:
                # defaults to all conditions if none were selected
                condition_names = gui.selected_condition_names or [ each.name for each in gui.conditions ]
            else:
                condition_names = [ each.name for each in gui.conditions ]
            
            return condition_names
            
        return wrapper
    
    def combined_wig_filtered_by_condition_input(panel, sizer, default_value=True):
        check_box_getter = create_check_box_getter(panel, sizer,
            label_text="Only Selected Conditions",
            default_value=default_value,
            tooltip_text="When checked, use the conditions table (on the left) to select which conditions to run this analysis on",
        )
        def wrapper(*args, **kwargs):
            is_checked = check_box_getter(*args, **kwargs)
            if is_checked:
                # defaults to all conditions if none were selected
                condition_names = gui.selected_condition_names or [ each.name for each in gui.conditions ]
            else:
                condition_names = [ each.name for each in gui.conditions ]
            
            return gui.combined_wigs[-1].with_only(condition_names=condition_names)
            
        return wrapper
    
    def combined_wig_filtered_by_sample_input(panel, sizer, default_value=True):
        check_box_getter = create_check_box_getter(panel, sizer,
            label_text="Only Selected Samples",
            default_value=default_value,
            tooltip_text="When checked, use the sample table (on the left) to select which samples to run this analysis on",
        )
        def wrapper(*args, **kwargs):
            is_checked = check_box_getter(*args, **kwargs)
            if is_checked:
                # defaults to all samples if none were selected
                selected_samples = gui.selected_samples or gui.samples
                wig_fingerprints = [ each.fingerprint for each in selected_samples ]
            else:
                wig_fingerprints = [ each.fingerprint for each in gui.samples ]
            
            return gui.combined_wigs[-1].with_only(wig_fingerprints=wig_fingerprints)
            
        return wrapper
    
    def create_run_button(panel, sizer, from_gui_function):
        from wx.lib.buttons import GenButton
        
        font = wx.Font(10, family = wx.FONTFAMILY_MODERN, style=0, weight=90, underline=False, faceName="", encoding=wx.FONTENCODING_DEFAULT)
  
        run_button = GenButton(panel, id=wx.ID_ANY, label="Run", style=wx.NO_BORDER, size=(-1, -1))
        run_button.SetBackgroundColour((78, 201, 176, 209))
        sizer.Add(10,20) # vertical padding
        sizer.Add(run_button, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL, border=gui_tools.default_padding)
        
        @gui_tools.bind_to(run_button, wx.EVT_BUTTON)
        def run(*args):
            # a workaround for python WX somehow clicking the add files button every time the run method is called
            def run_wrapper():
                with gui_tools.nice_error_log:
                    method_instance = from_gui_function(gui.frame)
                    gui.busy_running_method = True
                    if hasattr(method_instance, 'Run'):
                        method_instance.Run()
                gui.busy_running_method = False
                
            import threading
            thread = threading.Thread(target=run_wrapper())
            thread.setDaemon(True)
            thread.start()
                    
    def create_significance_choice_box(panel, sizer, default="HDI"):
        (
            signif_label,
            signif_wxobj,
            signif_sizer,
        ) = define_choice_box(
            panel,
            label_text="Significance method: ",
            options=["HDI","prob","BFDR","FWER"],
            tooltip_text="""Evaluate these various methods for determining significance of interactions.The options are:
            
            --HDI: significant genes are those for which the HDI does not overlap the ROPE
            --prob: significant genes are those with prob < 0.05, where ‘prob’ is probability that HDI overlaps the ROPE (default)
            --BFDR: significant genes are those with adjusted prob < 0.05, where prob is adjusted by the BFDR method
            --FWER: significant genes are those with adjusted prob < 0.05, where prob is adjusted by the FWER method""".replace("\n            ","\n"),
        )
        sizer.Add(signif_sizer, 1, wx.ALIGN_LEFT, gui_tools.default_padding)
        # return a value-getter
        signif_wxobj.SetSelection(signif_wxobj.FindString(default))
        return lambda *args: signif_wxobj.GetString(signif_wxobj.GetCurrentSelection())
    
    
