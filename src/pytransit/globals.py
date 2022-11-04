from random import random, seed
from os import getcwd, path
import sys

from pytransit.generic_tools import csv
from pytransit.generic_tools.lazy_dict import LazyDict, stringify, indent
from pytransit.generic_tools.named_list import named_list
from pytransit.generic_tools.misc import flatten_once, no_duplicates, singleton

debugging_enabled = True
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
        "Pre-Processing": {},
    })
    
    _annotation_path = "" if not debugging_enabled else f"{root_folder}/src/pytransit/data/genomes/H37Rv.prot_table"
    #_annotation_path = "" if not debugging_enabled else "H37Rv.prot_table"

    @property
    def annotation_path(self):
        from pytransit.specific_tools import transit_tools
        import os
        # validate it anytime the GUI tries to retrieve the annotation
        if not os.path.isfile(self._annotation_path):
            from pytransit.specific_tools import logging
            logging.error(f"Error: Annotation doesn't seem to be a file:{self._annotation_path}")
        return self._annotation_path
    
    @property
    def width(self):
        return self.frame.GetSize()[0]
    
    @property
    def height(self):
        return self.frame.GetSize()[1]
        
    @property
    def conditions(self):
        # only for latest combined wig
        return self.combined_wigs[-1].conditions
            
    @property
    def samples(self):
        # only for latest combined wig
        return self.combined_wigs[-1].samples
    
    @property
    def selected_samples(self): # this wrapper is here as an intentional design choice. Data access is done through a central place (this file) to allow for changing the implementation later without updating all the individual methods
        from pytransit.components.samples_area import get_selected_samples
        return get_selected_samples()
    
    @property
    def selected_condition_names(self): # this wrapper is here as an intentional design choice. Data access is done through a central place (this file) to allow for changing the implementation later without updating all the individual methods
        from pytransit.components.samples_area import get_selected_condition_names
        return get_selected_condition_names()
    
    @property
    def selected_conditions(self): # this wrapper is here as an intentional design choice. Data access is done through a central place (this file) to allow for changing the implementation later without updating all the individual methods
        from pytransit.components.samples_area import get_selected_condition_names
        names = get_selected_condition_names()
        return [ each for each in self.conditions if each.name in names ]
    
    @property
    def wigs_in_selected_conditions(self): # this wrapper is here as an intentional design choice. Data access is done through a central place (this file) to allow for changing the implementation later without updating all the individual methods
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
        from pytransit.components import samples_area
        return samples_area.add_wig_area_dropdown_option(*args, **kwargs)
    
    def add_condition_area_dropdown_option(self, *args, **kwargs):
        from pytransit.components import samples_area
        return samples_area.add_condition_area_dropdown_option(*args, **kwargs)
    
    def add_wig_area_button(self, *args, **kwargs):
        from pytransit.components import samples_area
        return samples_area.add_wig_area_button(*args, **kwargs)
    
    def add_result(self, *args, **kwargs):
        from pytransit.components import results_area
        return results_area.add(*args, **kwargs)
    
    debug_wx_python = True
    def debug_wx_if_needed(self):
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


if debugging_enabled:
    seed(0)
    import numpy
    numpy.random.seed(0)
