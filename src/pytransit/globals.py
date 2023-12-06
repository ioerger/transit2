from random import random, seed
from os import getcwd, path, getenv
import sys

import numpy
numpy.bool = bool

from pytransit.generic_tools import csv
from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.generic_tools.misc import flatten_once, no_duplicates, singleton

debugging_enabled = not not getenv("DEBUG")
root_folder       = path.join(path.dirname(__file__),"../../")

# 
# design
# 
    # @cli.add_command("resampling")
    # @gui.add_menu("Analysis")
    # @gui.add_wig_area_dropdown_option("Name")
    # @gui.add_condition_area_dropdown_option("Name")
    # @gui.add_results_button("Name")

# tools to make
    # a self-caching 'add file' button (remembers what files have been added)
    # a dynamically refreshing button area
    # a multi-select system

@singleton
class gui:
    is_active = len(sys.argv) <= 1
    frame = None
    busy_running_method = False

    combined_wigs = []
    
    menu_heirarchy = LazyDict({
        "File": {},
        "Pre-Processing": {},
    })
    

    @property
    def annotation_path(self):
        from pytransit.specific_tools import transit_tools
        import os
        _annotation_path = self.combined_wigs[-1].annotation_path
        # validate it anytime the GUI tries to retrieve the annotation
        if not os.path.isfile(_annotation_path):
            logging.error(f"Error: Annotation doesn't seem to be a file:{_annotation_path}")
        return _annotation_path
    
    @property
    def width(self):
        return self.frame.GetSize()[0]
    
    @property
    def height(self):
        return self.frame.GetSize()[1]
        
    @property
    def conditions(self):
        if len(self.combined_wigs):
            return self.combined_wigs[-1].conditions
        else:
            return []
            
    @property
    def samples(self):
        if len(self.combined_wigs):
            return self.combined_wigs[-1].samples
        else:
            return []
    
    @property
    def selected_samples(self): # this wrapper is here as an intentional design choice. Data access is done through a central place (this file) to allow for changing the implementation later without updating all the individual methods
        if not gui.is_active:
            return []
        from pytransit.components.samples_area import get_selected_samples
        return get_selected_samples()
    
    @property
    def selected_condition_names(self): # this wrapper is here as an intentional design choice. Data access is done through a central place (this file) to allow for changing the implementation later without updating all the individual methods
        if not gui.is_active:
            return []
        from pytransit.components.samples_area import get_selected_condition_names
        return get_selected_condition_names()
    
    @property
    def selected_conditions(self): # this wrapper is here as an intentional design choice. Data access is done through a central place (this file) to allow for changing the implementation later without updating all the individual methods
        if not gui.is_active:
            return []
        from pytransit.components.samples_area import get_selected_condition_names
        names = get_selected_condition_names()
        return [ each for each in self.conditions if each.name in names ]
    
    @property
    def wigs_in_selected_conditions(self): # this wrapper is here as an intentional design choice. Data access is done through a central place (this file) to allow for changing the implementation later without updating all the individual methods
        if not gui.is_active:
            return []
        selected_condition_names = set(self.selected_condition_names)
        selected_wigs = [ each for each in self.samples if len(set(each.condition_names) & selected_condition_names) > 0 ]
        return selected_wigs
    
    def add_menu(self, *args):
        """
        Example:
            @add_menu("Analysis", "Resampling")
            def on_menu_click():
                print("Menu item was clicked")
        """
        if not gui.is_active:
            return lambda func: func
        if len(args) <= 1:
            raise Exception(f'''When calling @gui.add_menu({args}), there needs to be at least two args (and all args need to be strings)''')
        
        parent_menu = gui.menu_heirarchy
        for each_name in args[:-1]:
            if not isinstance(parent_menu.get(each_name, None), dict):
                parent_menu[each_name] = LazyDict()
            # make sure all the child menus exist
            parent_menu = parent_menu[each_name]
        
        last_menu_name = args[-1]
        parent_menu[last_menu_name] = None
        
        def decorator(function_being_wrapped):
            from pytransit.specific_tools import gui_tools
            
            def wrapped_with_helper(*args, **kwargs):
                with gui_tools.nice_error_log:
                    return function_being_wrapped(*args, **kwargs)
            # attach the function
            parent_menu[last_menu_name] = wrapped_with_helper
            return function_being_wrapped
        return decorator

    def add_wig_area_dropdown_option(self, *args, **kwargs):
        if not gui.is_active:
            return lambda func: func
        from pytransit.components import samples_area
        return samples_area.add_wig_area_dropdown_option(*args, **kwargs)
    
    def add_condition_area_dropdown_option(self, *args, **kwargs):
        if not gui.is_active:
            return lambda func: func
        from pytransit.components import samples_area
        return samples_area.add_condition_area_dropdown_option(*args, **kwargs)
    
    def add_wig_area_button(self, *args, **kwargs):
        if not gui.is_active:
            return lambda func: func
        from pytransit.components import samples_area
        return samples_area.add_wig_area_button(*args, **kwargs)
    
    def add_result(self, *args, **kwargs):
        if not gui.is_active:
            return []
        from pytransit.components import results_area
        return results_area.add(*args, **kwargs)
    
    debug_wx_python = False
    def debug_wx_if_needed(self):
        if not gui.is_active:
            return []
        if self.debug_wx_python:
            import wx.lib.inspection
            wx.lib.inspection.InspectionTool().Show()

@singleton
class cli:
    # the order subcommands are shown can be defined here
    subcommands = {
        ("help",): "This should (and likely is) be replaced by a function elsewhere in the code"
    }
    
    def add_command(self, *args):
        """
        Example:
            @add_command("resampling")
            def on_run_command(args, kwargs):
                print(f"someone ran 'python ./src/transit resampling' with these args: {args} ")
        """
        def decorator(on_run_command_function):
            cli.subcommands[args] = on_run_command_function
            return on_run_command_function
        return decorator

@singleton
class logging:
    @staticmethod
    def log(message, *args, **kwargs):
        import inspect
        import os
        if gui.is_active:
            from pytransit.specific_tools.gui_tools import set_status
        else:
            def set_status(*args, **kwargs): pass
            
        
        message = f"{message} "+ " ".join([ f"{each}" for each in args])
        
        # get some context as to who is creating the message
        stack             = inspect.stack()
        caller_frame_info = stack[1]
        file_name         = ""
        caller_name       = ""
        try: file_name = os.path.basename(caller_frame_info.filename)
        except Exception as error: pass # sometimes the caller doesn't have a file name (ex: REPL)
        try: caller_name = caller_frame_info.function
        except Exception as error: pass # sometimes the caller doesn't have a function name (ex: lambda)
        
        # remove the .py extension
        if file_name[len(file_name)-3:len(file_name)] == ".py":
            file_name = file_name[0:len(file_name)-3]
        
        prefix = f"[{file_name}:{caller_name}()] "
        message = message.replace("\n", f"\n{prefix}")
        message = message.replace("\r", f"\r{prefix}")
        print(prefix+message, flush=True, **kwargs)
        if gui.is_active:
            set_status(message)
    
    @staticmethod
    def warn(*args, **kwargs):
        if gui.is_active:
            from pytransit.specific_tools.gui_tools import set_status
        else:
            def set_status(*args, **kwargs): pass
            
        import warnings
        string_output = _print_to_string(*args, **kwargs)
        warnings.warn("\n"+string_output, stacklevel=2)
        last_line = string_output.strip().split("\n")[-1]
        set_status(last_line)

    @staticmethod
    def error(*args, no_traceback=False, exit_code=1, **kwargs):
        import traceback
        if gui.is_active:
            from pytransit.specific_tools.gui_tools import set_status
        else:
            def set_status(*args, **kwargs): pass
        
        # this case is different from below because it uses stderr
        if no_traceback and not gui.is_active:
            import sys
            print(*args, file=sys.stderr, **kwargs)
            exit(exit_code)
        else:
            error_message = _print_to_string(*args, **kwargs)
            last_line = error_message.strip().split("\n")[-1]
            try:
                # this looks dumb but 'raise' is needed to get the traceback
                raise TransitError(error_message)
            except Exception as error:
                if not no_traceback:
                    traceback.print_exc()
                
                if gui.is_active:
                    set_status(last_line)
                    raise TransitError(error_message)
                else:
                    exit(exit_code)
    @staticmethod
    def progress_update(text, percent):
        if gui.is_active:
            from pytransit.components.parameter_panel import panel, progress_update
            progress_update(text, percent)
        else:
            print(text, end="\r")

class TransitError(Exception):
    pass
    
def _print_to_string(*args, **kwargs):
    from io import StringIO
    string_stream = StringIO()
    # dump to string
    print(*args, **{ "flush": True, **kwargs, "file":string_stream })
    output_str = string_stream.getvalue()
    string_stream.close()
    return output_str


if debugging_enabled:
    seed(0)
    import numpy
    numpy.random.seed(0)
