from pytransit.tools.transit_tools import wx, pub
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, stat_tools

default_label_size = (200, -1)
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
            self.wx_panel = wx.Panel(
                gui.frame,
                wx.ID_ANY,
                wx.DefaultPosition,
                wx.DefaultSize,
                #wx.Size(int(gui.frame.GetSize()[0]/2), wx.DefaultSize[1]),
                wx.TAB_TRAVERSAL,
            )
            self.main_sizer = wx.BoxSizer(wx.VERTICAL)
            return (self.wx_panel, self.main_sizer)
        
        def __exit__(self, _, error, traceback_obj):
            if error is not None:
                import traceback
                print(''.join(traceback.format_tb(traceback_obj)))
                frame = gui.frame
                if frame and hasattr(frame, "status_bar"):
                    frame.status_bar.SetStatusText("Error: "+str(error.args))
            else:
                from pytransit.components import parameter_panel
                parameter_panel.set_panel(self.wx_panel)
                self.wx_panel.SetSizer(self.main_sizer)
                self.wx_panel.Layout()
                self.main_sizer.Fit(self.wx_panel)
                gui.frame.Layout()
    
    def create_button(panel, sizer, *, label):
        run_button = wx.Button(
            panel,
            wx.ID_ANY,
            label,
            wx.DefaultPosition,
            wx.DefaultSize,
            0,
        )
        sizer.Add(run_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        def decorator(func):
            @gui_tools.bind_to(run_button, wx.EVT_BUTTON)
            def wrapper(*args,**kwargs):
                with gui_tools.nice_error_log:
                    return func(*args,**kwargs)
            return wrapper
        return decorator

    

    def create_default_pathway_button(panel, sizer, *, button_label, tooltip_text="", popup_title=""):
        import csv
        COG_orgs = []
        with open(root_folder+"src/pytransit/data/cog-20.org.csv") as file_obj:
            reader_obj = csv.reader(file_obj)
            for row in reader_obj:
                COG_orgs.append(row[1])

        """
            Example:
                file_path_getter = create_file_input(self.panel, main_sizer, button_label="Add context file", allowed_extensions='All files (*.*)|*.*')
                file_path_or_none = file_path_getter()
        """
        from os.path import basename
        row_sizer = wx.BoxSizer(wx.HORIZONTAL)
        if True:
            # 
            # tooltip
            # 
            if tooltip_text:
                from pytransit.components.icon import InfoIcon
                row_sizer.Add(
                    InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
                    0,
                    wx.ALIGN_CENTER_VERTICAL,
                    gui_tools.default_padding,
                )
            # 
            # button
            # 
            if True:
                popup_button = wx.Button(
                    panel,
                    wx.ID_ANY,
                    button_label,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                # whenever the button is clicked, popup
                organism_pathway_text = None
                organism_pathway = None
                @gui_tools.bind_to(popup_button, wx.EVT_BUTTON)
                def when_button_clicked(*args,**kwargs):
                    nonlocal organism_pathway
                    win = wx.Dialog(panel,wx.FRAME_FLOAT_ON_PARENT)
                    popup_sizer = wx.BoxSizer(wx.VERTICAL)
                    win.SetSizer(popup_sizer)

                    pathway_label_text= wx.StaticText(win, wx.ID_ANY, label="Select A Pathway Type : ", style=wx.ALIGN_LEFT)
                    popup_sizer.Add(pathway_label_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                    pathway_type = wx.ComboBox(win,choices = ["Sanger", "COG", "COG_20" ,"GO", "KEGG"])
                    popup_sizer.Add(pathway_type,wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)

                    select_btn = wx.Button(win, wx.ID_OK, label = "Select", size = (50,20), pos = (75,50))
                    popup_sizer.Add(select_btn,wx.EXPAND, gui_tools.default_padding)

                    win.Layout()
                    popup_sizer.Fit(win)
                    selected_path = win.ShowModal()

                    if selected_path == wx.ID_OK:
                        pathway_type_selected = pathway_type.GetValue()

                        if pathway_type_selected== "COG_20":
                            organism_label_text= wx.StaticText(win, wx.ID_ANY, label="Select An Organism : ", style=wx.ALIGN_LEFT)
                            popup_sizer.Add(organism_label_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                            organism = wx.ComboBox(win,choices = COG_orgs)
                            popup_sizer.Add(organism,wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                        else:
                            organism_label_text= wx.StaticText(win, wx.ID_ANY, label="Select An Organism : ", style=wx.ALIGN_LEFT)
                            popup_sizer.Add(organism_label_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                            organism = wx.ComboBox(win,choices = ["H37Rv", "Smeg"])
                            popup_sizer.Add(organism,wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)

                        ok_btn = wx.Button(win, wx.ID_OK, label = "Ok", size = (50,20), pos = (75,50))
                        popup_sizer.Add(ok_btn,wx.EXPAND, gui_tools.default_padding)

                        win.Layout()
                        popup_sizer.Fit(win)
                        res = win.ShowModal()
                        if res == wx.ID_OK:
                            organism_pathway = "-".join([organism.GetValue(),pathway_type_selected])
                            organism_pathway_text.SetLabel(basename(organism_pathway or ""))
                        win.Destroy()
                    

            row_sizer.Add(popup_button, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)

            organism_pathway_text= wx.StaticText(panel, wx.ID_ANY, label="", style=wx.ALIGN_LEFT)
            row_sizer.Add(organism_pathway_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
        
        sizer.Add(row_sizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        return lambda *args, **kwargs: organism_pathway



    def create_file_input(panel, sizer, *, button_label, tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*'):
        """
            Example:
                file_path_getter = create_file_input(self.panel, main_sizer, button_label="Add context file", allowed_extensions='All files (*.*)|*.*')
                file_path_or_none = file_path_getter()
        """
        from os.path import basename
        row_sizer = wx.BoxSizer(wx.HORIZONTAL)
        if True:
            # 
            # tooltip
            # 
            if tooltip_text:
                from pytransit.components.icon import InfoIcon
                row_sizer.Add(
                    InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
                    0,
                    wx.ALIGN_CENTER_VERTICAL,
                    gui_tools.default_padding,
                )
            
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
            row_sizer.Add(add_file_button, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
            
            # 
            # Text
            # 
            file_text = wx.StaticText(panel, wx.ID_ANY, label="", style=wx.ALIGN_LEFT)
            row_sizer.Add(file_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
        
        sizer.Add(row_sizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        return lambda *args, **kwargs: the_file_path
    
    def create_persistent_file_input(panel, sizer, *, name, button_label, tooltip_text="", popup_title="", default_folder=None, default_file_name="", allowed_extensions='All files (*.*)|*.*'):
        """
            Example:
                file_path_getter = create_persistent,file_input(self.panel, main_sizer, button_label="Add context file", allowed_extensions='All files (*.*)|*.*')
                file_path_or_none = file_path_getter()
        """
        ##todo : HANDLE NAME
        from os.path import basename
        row_sizer = wx.BoxSizer(wx.HORIZONTAL)
        if True:
            # 
            # tooltip
            # 
            if tooltip_text:
                from pytransit.methods.analysis_base import InfoIcon
                row_sizer.Add(
                    InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
                    0,
                    wx.ALIGN_CENTER_VERTICAL,
                    gui_tools.default_padding,
                )
            
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
            row_sizer.Add(add_file_button, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
            
            # 
            # Text
            # 
            file_text = wx.StaticText(panel, wx.ID_ANY, label="", style=wx.ALIGN_LEFT)
            row_sizer.Add(file_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
        
        sizer.Add(row_sizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
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
            widget_size = (100, -1)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(panel, wx.ID_ANY, label_text, wx.DefaultPosition, label_size, 0)
        label.Wrap(-1)
        choice_box = wx.Choice(panel, wx.ID_ANY, wx.DefaultPosition, widget_size, options, 0 )
        choice_box.SetSelection(0)
        sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
        sizer.Add(choice_box, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
        sizer.Add(
            InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
            0,
            wx.ALIGN_CENTER_VERTICAL,
            gui_tools.default_padding,
        )
        return (label, choice_box, sizer)
    
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
        
        sizer.Add(inner_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)

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
        sizer.Add(inner_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
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

            sizer = wx.BoxSizer(wx.HORIZONTAL)
            label = wx.StaticText(panel, wx.ID_ANY, label_text, wx.DefaultPosition, label_size, 0)
            label.Wrap(-1)
            text_box = wx.TextCtrl(panel, wx.ID_ANY, f"{default_value}", wx.DefaultPosition, widget_size, 0)
            
            sizer.Add(label, 0,  wx.ALL|wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
            sizer.Add(text_box, 0,  wx.ALL|wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
            sizer.Add(InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text), 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
            sizer.Layout()
            
            
            return label, text_box, sizer
            
    def create_check_box_getter(panel, sizer, label_text="", default_value=False, tooltip_text="", widget_size=None):
        from pytransit.components.icon import InfoIcon
        if not widget_size:
            widget_size = (-1, -1)
        
        inner_sizer = wx.BoxSizer(wx.HORIZONTAL)
        check_box   = wx.CheckBox(panel, label=label_text, size=widget_size)
        
        check_box.SetValue(default_value)
        
        inner_sizer.Add(check_box, 0, wx.ALIGN_CENTER_VERTICAL, gui_tools.default_padding)
        inner_sizer.Add(
            InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
            0,
            wx.ALIGN_CENTER_VERTICAL,
            gui_tools.default_padding,
        )
        
        sizer.Add(inner_sizer, 1, wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        
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
        sizer.Add(wrapper_sizer, 1, wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
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

# 
# 
# specific
# 
# 
if True:
    def create_preview_loess_button(panel, sizer, conditions_getter=None, wig_ids_getter=None):
        conditions_getter = conditions_getter or (lambda *args, **kwargs: None)
        wig_ids_getter    = wig_ids_getter    or (lambda *args, **kwargs: None)
        @create_button(panel, sizer, label="Preview LOESS fit")
        def when_loess_preview_clicked(event):
            with gui_tools.nice_error_log:
                import numpy
                import matplotlib
                import matplotlib.pyplot as plt
                from pytransit.tools import stat_tools
                from pytransit.globals import gui, cli, root_folder, debugging_enabled
                from pytransit.tools.tnseq_tools import Wig
                
                # 
                # determine selection method
                # 
                conditions = conditions_getter()
                wig_ids    = wig_ids_getter()
                use_selected = False
                if conditions is None or len(conditions) == 0 and (wig_ids is None or len(wig_ids) == 0):
                    use_selected = True
                    if not gui.selected_samples:
                        # NOTE: was a popup
                        logging.error("Need to select at least one control or experimental dataset.")
                
                wig_objects = gui.samples
                
                #
                # get read_counts and positions
                # 
                read_counts_per_wig, position_per_line = transit_tools.gather_sample_data_for(conditions=conditions, wig_ids=wig_ids, selected_samples=use_selected)
                number_of_wigs, number_of_lines = read_counts_per_wig.shape # => number_of_lines = len(position_per_line)
                window = 100
                for each_path_index in range(number_of_wigs):

                    number_of_windows = int(number_of_lines / window) + 1  # python3 requires explicit rounding to int
                    x_w = numpy.zeros(number_of_windows)
                    y_w = numpy.zeros(number_of_windows)
                    for window_index in range(number_of_windows):
                        x_w[window_index] = window * window_index
                        y_w[window_index] = sum(read_counts_per_wig[each_path_index][window * window_index : window * (window_index + 1)])
                    
                    y_smooth = stat_tools.loess(x_w, y_w, h=10000)
                    plt.plot(x_w, y_w, "g+")
                    plt.plot(x_w, y_smooth, "b-")
                    plt.xlabel("Genomic Position (TA sites)")
                    plt.ylabel("Reads per 100 insertion sites")
                    
                    plt.title("LOESS Fit - %s" % wig_objects[each_path_index].id)
                    plt.show()
    
    def create_normalization_input(panel, sizer, default="TTR"):
        (
            label,
            normalization_wxobj,
            normalization_choice_sizer,
        ) = define_choice_box(
            panel,
            label_text="Normalization: ",
            options=[
                "TTR",
                "nzmean",
                "totreads",
                "zinfnb",
                "quantile",
                "betageom",
                "nonorm",
            ],
            tooltip_text="Choice of normalization method. The default choice, 'TTR', normalizes datasets to have the same expected count (while not being sensative to outliers). Read documentation for a description other methods. ",
        )
        sizer.Add(normalization_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        # return a value-getter
        normalization_wxobj.SetSelection(normalization_wxobj.FindString(default))
        return lambda *args: normalization_wxobj.GetString(normalization_wxobj.GetCurrentSelection())
    
    def create_condition_choice(panel, sizer, *, label_text, tooltip_text="choose condition"):
        (
            label,
            ref_condition_wxobj,
            ref_condition_choice_sizer,
        ) = define_choice_box(
            panel,
            label_text=label_text,
            options=[x.name for x in gui.conditions],
            tooltip_text=tooltip_text,
        )
        sizer.Add(ref_condition_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        return lambda *args: ref_condition_wxobj.GetString(ref_condition_wxobj.GetCurrentSelection())
    
    def create_reference_condition_input(panel, sizer):
        (
            label,
            ref_condition_wxobj,
            ref_condition_choice_sizer,
        ) = define_choice_box(
            panel,
            label_text="Ref Condition:",
            options=[ "[None]" ] + [ each.name for each in gui.conditions ],
            tooltip_text="which condition(s) to use as a reference for calculating LFCs (comma-separated if multiple conditions)",
        )
        sizer.Add(ref_condition_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        return lambda *args: ref_condition_wxobj.GetString(ref_condition_wxobj.GetCurrentSelection())
    
    def create_condition_input(panel, sizer, label_text="Condition", tooltip_text="choose condition"):
        (
            label,
            ref_condition_wxobj,
            ref_condition_choice_sizer,
        ) = define_choice_box(
            panel,
            label_text=label_text,
            options=[ "[None]" ] + [x.name for x in gui.conditions],
            tooltip_text=tooltip_text,
        )
        sizer.Add(ref_condition_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        return lambda *args: ref_condition_wxobj.GetString(ref_condition_wxobj.GetCurrentSelection())
    
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
        sizer.Add(wrapper_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        
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
        sizer.Add(wrapper_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        
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
            tooltip_text="which condition(s) to use as the control group",
            label_size=(200, 20),
        )
        sizer.Add(ref_condition_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
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
            tooltip_text="which condition(s) to use as the experimental group",
            label_size=(200, 20),
        )
        sizer.Add(ref_condition_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        return lambda *args: ref_condition_wxobj.GetString(ref_condition_wxobj.GetCurrentSelection())
    
    def create_n_terminus_input(panel, sizer):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text="Ignore N-Terminus %:",
            default_value="0",
            tooltip_text="Ignores a fraction of the ORF, beginning at the N-terminal end. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",    
        )
        return lambda *args: float(get_text())
    
    def create_c_terminus_input(panel, sizer):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text="Ignore C-Terminus %:",
            default_value="0",
            tooltip_text="Ignores a fraction of the ORF, beginning at the C-terminal end. Useful for ignoring read-counts that may occur at the terminal ends, even though they do not truly disrupt a genes function.",    
        )
        return lambda *args: float(get_text())
    
    def create_pseudocount_input(panel, sizer, default_value="5"):
        # 
        # text input: Pseudocount
        # 
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text="Pseudocount:",
            default_value=default_value,
            tooltip_text="Pseudo-counts used in calculating log-fold-change. Useful to dampen the effects of small counts which may lead to deceptively high LFC.",    
        )
        return lambda *args: float(get_text())
    
    def create_pseudocount_input(panel, sizer, default_value="5"):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text="Pseudocount:",
            default_value=default_value,
            tooltip_text="Pseudo-counts used in calculating log-fold-change. Useful to dampen the effects of small counts which may lead to deceptively high LFC.",    
        )
        return lambda *args: int(get_text())
    
    def create_alpha_input(panel, sizer, default_value="1000"):
        get_text = create_text_box_getter(
            panel,
            sizer,
            label_text="Alpha:",
            default_value=default_value,
            tooltip_text="FIXME: needs summary",    
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
    
    def create_run_button(panel, sizer, from_gui_function):
        run_button = wx.Button(
            panel,
            wx.ID_ANY,
            "Run",
            wx.DefaultPosition,
            wx.DefaultSize,
            0,
        )
        sizer.Add(run_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        
        
        @gui_tools.bind_to(run_button, wx.EVT_BUTTON)
        def run(*args):
            # a workaround for python WX somehow clicking the add files button every time the run method is called
            def run_wrapper():
                gui.busy_running_method = True
                try:
                    method_instance.Run()
                except Exception as error:
                    pass
                gui.busy_running_method = False
                
            import threading
            with gui_tools.nice_error_log:
                method_instance = from_gui_function(gui.frame)
                if method_instance:
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
            tooltip_text="tooltip",  # FIXME: fill in explanation...
        )
        sizer.Add(signif_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        # return a value-getter
        signif_wxobj.SetSelection(signif_wxobj.FindString(default))
        return lambda *args: signif_wxobj.GetString(signif_wxobj.GetCurrentSelection())
