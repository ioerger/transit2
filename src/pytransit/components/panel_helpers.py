from pytransit.transit_tools import wx, pub
from pytransit.core_data import universal

LABEL_SIZE = (100, -1)
WIDGET_SIZE = (100, -1)
default_padding = 5


# 
# 
# generic
# 
# 
if True:
    def define_choice_box(
        panel,
        label_text="",
        widget_choice=[""],
        tooltip_text="",
        lab_size=None,
        widget_size=None,
    ):
        from pytransit.analysis.base import InfoIcon
        
        if not lab_size:
            lab_size = (100, -1)
        if not widget_size:
            widget_size = (100, -1)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(panel, wx.ID_ANY, label_text, wx.DefaultPosition, lab_size, 0)
        label.Wrap(-1)
        choice_box = wx.Choice(panel, wx.ID_ANY, wx.DefaultPosition, widget_size, widget_choice, 0 )
        choice_box.SetSelection(0)
        sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(choice_box, 0, wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(
            InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
            0,
            wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        return (label, choice_box, sizer)


    def define_text_box(
            panel,
            label_text="",
            widget_text="",
            tooltip_text="",
            lab_size=None,
            widget_size=None,
        ):
            from pytransit.analysis.base import InfoIcon
            if not lab_size:
                lab_size = LABEL_SIZE
            if not widget_size:
                widget_size = WIDGET_SIZE

            sizer = wx.BoxSizer(wx.HORIZONTAL)
            label = wx.StaticText(
                panel, wx.ID_ANY, label_text, wx.DefaultPosition, lab_size, 0
            )
            label.Wrap(-1)
            text_box = wx.TextCtrl(
                panel, wx.ID_ANY, widget_text, wx.DefaultPosition, widget_size, 0
            )
            sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL, 5)
            sizer.Add(text_box, 0, wx.ALIGN_CENTER_VERTICAL, 5)
            sizer.Add(
                InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
                0,
                wx.ALIGN_CENTER_VERTICAL,
                5,
            )
            return (label, text_box, sizer)

# 
# 
# specific
# 
# 
if True:
    def create_normalization_dropdown(panel, sizer):
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
        sizer.Add(normalization_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, default_padding)
        # return a value-getter
        return lambda *args: normalization_wxobj.GetString(normalization_wxobj.GetCurrentSelection())
    
    def create_reference_condition_dropdown(panel, sizer):
        (
            label,
            ref_condition_wxobj,
            ref_condition_choice_sizer,
        ) = define_choice_box(
            panel,
            "Ref Condition:",
            universal.session_data.condition_names,
            "which condition(s) to use as a reference for calculating LFCs (comma-separated if multiple conditions)",
        )
        sizer.Add(ref_condition_choice_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, default_padding)
        return lambda *args: ref_condition_wxobj.GetString(ref_condition_wxobj.GetCurrentSelection())
    
    def create_include_condition_list(panel, sizer):
        (
            _,
            wxobj,
            wrapper_sizer,
        ) = define_text_box(
            panel=panel,
            label_text="Include\nConditions\n",
            widget_text="",
            tooltip_text="comma seperated list (default=all)",
        )
        sizer.Add(wrapper_sizer, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, default_padding)
        
        def get_value(*args):
            as_list = ",".split(wxobj.GetValue())
            if len(as_list) == 0:
                return universal.session_data.condition_names,
        
        return get_value