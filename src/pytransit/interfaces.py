from random import random, seed
from os import getcwd, path
import sys

import pytransit.basics.csv as csv
from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.basics.misc import flatten_once, no_duplicates, singleton

debugging_enabled = True
root_folder       = path.join(path.dirname(__file__),"../../")
if debugging_enabled:
    seed(0)
    
# # new design

# gui # instead of universal
# @gui.add_menu("Analysis")
# @gui.add_sample_button("Name")

# @cli.from_args("name")


@singleton
class gui:
    is_active = len(sys.argv) <= 1
    frame = None
    busy_running_method = False

    annotation_path = "" if not debugging_enabled else f"{getcwd()}/src/pytransit/genomes/H37Rv_dev.prot_table"
    combined_wigs = []
    
    menu_heirarchy = LazyDict()
    
    @property
    def conditions(self):
        return no_duplicates(flatten_once(each_combined_wig.conditions for each_combined_wig in self.combined_wigs))
            
    @property
    def samples(self):
        return no_duplicates(flatten_once(each_combined_wig.samples for each_combined_wig in self.combined_wigs))
    
    @property
    def selected_samples(self):
        from pytransit.components.samples_area import get_selected_samples
        return get_selected_samples()
    
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
            from pytransit.tools import gui_tools
            
            def wrapped_with_helper(*args, **kwargs):
                with gui_tools.nice_error_log:
                    return function_being_wrapped(*args, **kwargs)
            # attach the function
            parent_menu[last_menu_name] = wrapped_with_helper
            return function_being_wrapped
        return decorator
