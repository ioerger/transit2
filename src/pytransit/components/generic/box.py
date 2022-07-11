from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools

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
    def __init__(self, orientation=None, min_size=None, max_size=None, children=None):
        window       = gui_tools.window
        children     = children or []
        orientation  = orientation if orientation is not None else options_for.orientation.vertical
        
        # 
        # wx_object
        # 
        wx_object = wx.BoxSizer(orient=orientation)
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
            self.add(children)
        
    def add(self, component, proportion=1, side=None, expand=None, horizontal_alignment=None):
        side                 = 0b00000000                            if side                 is None else options_for.sides[side]
        expand               = 0b00000000                            if expand               is None else options_for.expand[expand]
        horizontal_alignment = options_for.horizontal_alignment.left if horizontal_alignment is None else options_for.horizontal_alignment[horizontal_alignment]
        
        self._state.children.append(component)
        
        wx_object = component.wx_object if hasattr(component, "wx_object") else component
        
        self.wx_object.Add(
            wx_object,
            proportion,
            flag=side | expand | horizontal_alignment,
            border=0,
        )
        return self
        
    @property
    def length(self):
        return len(self.children)
    
    # TODO: remove
    # TODO: insert
    
    def __len__(self):
        return self.length
    
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)

class Column:
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, min_size=None, max_size=None, children=None):
        super(Column, self).__init__(orientation="vertical", min_size=min_size, max_size=max_size, children=children)

class Row:
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, min_size=None, max_size=None, children=None):
        super(Row, self).__init__(orientation="horizontal", min_size=min_size, max_size=max_size, children=children)