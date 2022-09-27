from random import random
from os import getcwd, path

import pytransit.basics.csv as csv
from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.basics.misc import flatten_once, no_duplicates, singleton

# TODO:
    # should probably split this into a universal and a gui_tools.global
    # because the from_args isn't going to use a lot of these things
    # although .Run() might use some

# 
# universal
# 
    # tools and methods ideally would use this interface instead of accessing individual GUI components
    # that way, the GUI can change as much as needed, so long as it maintains the universal data interface

universal = LazyDict(
    interface=None, # "gui" or "console"
    frame=None,
    debugging_enabled=False,
    session_data=None,
    busy_running_method=False,
    root_folder = path.join(path.dirname(__file__),"../../")
)

@singleton
class SessionData(LazyDict):
    annotation_path = "" if not universal.debugging_enabled else f"{getcwd()}/src/pytransit/genomes/H37Rv_dev.prot_table"
    combined_wigs = []
    
    @property
    def conditions(self):
        return no_duplicates(flatten_once(each_combined_wig.conditions for each_combined_wig in self.combined_wigs))
            
    @property
    def samples(self):
        return no_duplicates(flatten_once(each_combined_wig.samples for each_combined_wig in self.combined_wigs))
    
    @property
    def selected_samples(self):
        if universal.interface == "gui":
            from pytransit.components.samples_area import get_selected_samples
            return get_selected_samples()
        else:
            # currently (Aug 2022) CLI doesn't use or designate a form of selected samples
            return self.samples

universal.session_data = SessionData