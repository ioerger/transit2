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
    def __init__(self, orientation=None, background_color=None, min_size=None, max_size=None, children=None):
        window       = gui_tools.window
        children     = children or []
        orientation  = options_for.orientation[orientation] if orientation is not None else options_for.orientation.vertical
        
        # 
        # wx_object
        # 
        wx_object = wx.BoxSizer(orient=orientation)
        panel = wx.Panel(window)
        if background_color: panel.SetBackgroundColour(background_color)
        sizer = wx.BoxSizer(orient=orientation)
        if max_size: sizer.SetMaxSize(wx.Size(*max_size))
        if min_size: sizer.SetMinSize(wx.Size(*min_size))
        wx_object.Add(panel)
        
        self.wx_object = wx_object
        self.events = LazyDict(
            # None
        )
        self._state = LazyDict(
            panel=panel,
            sizer=sizer,
            children=[],
        )
        for each in children:
            self.add(each)
        
        self.refresh()
        
    def add(self, component, proportion=1, side=None, expand=None, horizontal_alignment=None):
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
            
        self._state.sizer.Add(
            wx_object,
            proportion,
            flag=side | expand | horizontal_alignment,
            border=0,
        )
        return self
        
    @property
    def length(self):
        return len(self.children)
    
    def refresh(self):
        self._state.panel.SetSizer(self._state.sizer)
        self._state.panel.Layout()
        self._state.panel.Center(wx.BOTH)
    
    # TODO: remove
    # TODO: insert
    
    def __len__(self):
        return self.length
    
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        self.refresh()
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)

class Column(Box):
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, background_color=None, min_size=None, max_size=None, children=None):
        super(Column, self).__init__(orientation="vertical", background_color=background_color, min_size=min_size, max_size=max_size, children=children)

class Row(Box):
    """
        Overview:
            self.wx_object
            self.add(component)
    """
    def __init__(self, background_color=None, min_size=None, max_size=None, children=None):
        super(Row, self).__init__(orientation="horizontal", background_color=background_color, min_size=min_size, max_size=max_size, children=children)