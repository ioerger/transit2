from pytransit.transit_tools import wx, pub

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
    choiceBox = wx.Choice(panel, wx.ID_ANY, wx.DefaultPosition, widget_size, widget_choice, 0 )
    choiceBox.SetSelection(0)
    sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL, 5)
    sizer.Add(choiceBox, 0, wx.ALIGN_CENTER_VERTICAL, 5)
    sizer.Add(
        InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
        0,
        wx.ALIGN_CENTER_VERTICAL,
        5,
    )
    return (label, choiceBox, sizer)

def create_normalization_dropdown(panel):
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
    return normalization_choice_sizer, lambda : normalization_wxobj.GetString(normalization_wxobj.GetCurrentSelection())