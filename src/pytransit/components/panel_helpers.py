from pytransit.tools.transit_tools import wx, pub
from pytransit.universal_data import universal
import pytransit.tools.gui_tools as gui_tools

default_label_size = (-1, -1)
default_widget_size = (100, -1)

# 
# 
# generic
# 
# 
if True:
    def make_panel():
        return wx.Panel(
            universal.frame,
            wx.ID_ANY,
            wx.DefaultPosition,
            wx.DefaultSize,
            wx.TAB_TRAVERSAL,
        )
    
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
    
    def define_choice_box(
        panel,
        label_text="",
        options=[""],
        tooltip_text="",
        label_size=None,
        widget_size=None,
    ):
        from pytransit.methods.analysis_base import InfoIcon
        
        if not label_size:
            label_size = (100, -1)
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


    def define_text_box(
            panel,
            label_text="",
            default_value="",
            tooltip_text="",
            label_size=None,
            widget_size=None,
        ):
            from pytransit.methods.analysis_base import InfoIcon
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
        from pytransit.methods.analysis_base import InfoIcon
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
            default_value=default_value,
            tooltip_text=tooltip_text,
        )
        sizer.Add(wrapper_sizer, 1, wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        return lambda *args: wxobj.GetValue()

# 
# 
# specific
# 
# 
if True:
    def create_normalization_input(panel, sizer, default="TTR"):
        (
            label,
            normalization_wxobj,
            normalization_choice_sizer,
        ) = define_choice_box(
            panel,
            "Normalization: ",
            [
                "TTR",
                "nzmean",
                "totreads",
                "zinfnb",
                "quantile",
                "betageom",
                "nonorm",
            ],
            "Choice of normalization method. The default choice, 'TTR', normalizes datasets to have the same expected count (while not being sensative to outliers). Read documentation for a description other methods. ",
        )
        sizer.Add(normalization_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        # return a value-getter
        normalization_wxobj.SetSelection(normalization_wxobj.FindString(default))
        return lambda *args: normalization_wxobj.GetString(normalization_wxobj.GetCurrentSelection())
    
    def create_reference_condition_input(panel, sizer):
        (
            label,
            ref_condition_wxobj,
            ref_condition_choice_sizer,
        ) = define_choice_box(
            panel,
            "Ref Condition:",
            [ "[None]" ] + [ each.name for each in universal.session_data.conditions ],
            "which condition(s) to use as a reference for calculating LFCs (comma-separated if multiple conditions)",
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
                return [ each.name for each in universal.session_data.conditions ]
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
            "Control Condition:",
            [ "[None]" ] + [ each.name for each in universal.session_data.conditions ],
            "which condition(s) to use as the control group",
            label_size=(150, 20),
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
            "Experimental Condition:",
            [ "[None]" ] + [ each.name for each in universal.session_data.conditions ],
            "which condition(s) to use as the experimental group",
            label_size=(160, 20),
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
    
    def create_run_button(panel, sizer):
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
                universal.busy_running_method = True
                try:
                    method_instance.Run()
                except Exception as error:
                    pass
                universal.busy_running_method = False
                
            import threading
            with gui_tools.nice_error_log:
                method_instance = universal.selected_method.method.from_gui(universal.frame)
                if method_instance:
                    thread = threading.Thread(target=run_wrapper())
                    thread.setDaemon(True)
                    thread.start()
