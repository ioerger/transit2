from random import random
from os import getcwd

import pytransit.basics.csv as csv
from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.basics.misc import flatten_once, no_duplicates


class SessionData(LazyDict):
    annotation_path = f"{getcwd()}/src/pytransit/genomes/H37Rv_dev.prot_table" # FIXME: default value is for debugging only
    standalone_wigs = []
    combined_wigs = []
    
    @property
    def conditions(self): return no_duplicates(flatten_once(each.conditions for each in self.combined_wigs))
    @property
    def samples(self): return flatten_once(each.samples for each in self.combined_wigs)
    @property
    def selected_samples(self):
        from pytransit.components.samples_area import sample_table
        return [
            each["__wig_obj"]
                for each in sample_table.selected_rows
        ]

universal = LazyDict(
    frame=None,
    interface=None, # "gui" or "console"
    session_data=SessionData(),
    busy_running_method=False,
)