from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
from pytransit.specific_tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, basename
from pytransit.specific_tools import gui_tools

if HAS_WX:
    options_for = LazyDict(
        orientation = LazyDict(
            vertical=wx.VERTICAL,
            horizontal=wx.HORIZONTAL,
        ),
        sides = LazyDict(
            top = wx.TOP,
            bottom = wx.BOTTOM,
            left = wx.LEFT,
            right = wx.RIGHT,
            all = wx.ALL,
        ),
        expand = LazyDict({
            True: wx.EXPAND,
            "shaped": wx.SHAPED,
        }),
        vertical_alignment = LazyDict(
            top=wx.ALIGN_TOP,
            bottom=wx.ALIGN_BOTTOM,
            center=wx.ALIGN_CENTER_VERTICAL,
        ),
        horizontal_alignment = LazyDict(
            left=wx.ALIGN_LEFT,
            right=wx.ALIGN_RIGHT,
            center=wx.ALIGN_CENTER_HORIZONTAL,
        ),
    )

class Box:
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, orientation=None, background_color=None, min_size=None, max_size=None, children=None):
        children     = children or []
        orientation  = options_for.orientation[orientation] if orientation is not None else options_for.orientation.vertical
        
        # 
        # wx_object
        # 
        wx_object = wx.BoxSizer(orient=orientation)
        if background_color: wx_object.SetBackgroundColour(gui_tools.color.green)
        if max_size: wx_object.SetMaxSize(wx.Size(*max_size))
        if min_size: wx_object.SetMinSize(wx.Size(*min_size))
        
        self.wx_object = wx_object
        self.events = LazyDict(
            # None
        )
        self._state = LazyDict(
            children=[],
        )
        for each in children:
            self.add(each)
        
    def add(self, component, proportion=1, side=None, expand=None, horizontal_alignment=None, border=0):
        # handle multiple
        if isinstance(component, (list, tuple)):
            for each in component:
                self.add(each, proportion, side, expand, horizontal_alignment)
            return self
        
        side                 = 0b00000000                            if side                 is None else options_for.sides[side]
        expand               = 0b00000000                            if expand               is None else options_for.expand[expand]
        horizontal_alignment = options_for.horizontal_alignment.left if horizontal_alignment is None else options_for.horizontal_alignment[horizontal_alignment]
        
        self._state.children.append(component)
        
        wx_object = component.wx_object if hasattr(component, "wx_object") else component
        
        # some children need a parent before they can be created...for some reason
        if callable(wx_object):
            wx_object = wx_object(self.wx_object)
        
        # help(self.wx_object.Add)    
        self.wx_object.Add(
            wx_object,
            int(proportion),
            flag=(side or 0) | (expand or 0) | (horizontal_alignment or 0),
            border=border,
        )
        return self
        
    @property
    def length(self):
        return len(self.children)
    
    def __len__(self):
        return self.length
    
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)

class Column(Box):
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, min_size=None, max_size=None, children=None):
        super(Column, self).__init__(orientation="vertical", min_size=min_size, max_size=max_size, children=children)

class Row(Box):
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, min_size=None, max_size=None, children=None):
        super(Row, self).__init__(orientation="horizontal", min_size=min_size, max_size=max_size, children=children)