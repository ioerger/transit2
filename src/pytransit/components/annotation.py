from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
import pytransit.gui_tools as gui_tools

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button

class Annotation:
    """
        Overview:
            self.wx_object
            self.add(python_obj)
    """
    def __init__(self):
        window = gui_tools.window
        
        # 
        # data
        # 
        self.annotation_file = None
        
        # 
        # define component
        # 
        prefix         = Text("Annotation File:", bold=True)
        file_path_text = Text("[no file selected]")
        button = Button(
            text="Select File",
        )
        self._as_row = Row(children=[
            prefix,
            file_path_text,
            button,
        ])
        
        # 
        # define callbacks
        # 
        @button.events.on_click
        def when_clicked(event):
            files = gui_tools.ask_for_files("Select Annotation file (.prot_table or .gff3)")
            if files and len(files) == 1:
                self.annotation_file = files[0]
                # update the UI with the file content
                file_path_text.content = basename(self.annotation_file)
                # call all the callbacks
                for each in self._state.on_file_callbacks:
                    each(event)
                
        
        # 
        # standardize
        # 
        self.wx_object = self._as_row.wx_object
        self.events = LazyDict(
            on_file=lambda func: on_file_callbacks.append(func),
        )
        self._state = LazyDict(
            on_file_callbacks=[],
        )
    
    def __enter__(self):
        return self
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            gui_tools.handle_traceback(traceback_obj)