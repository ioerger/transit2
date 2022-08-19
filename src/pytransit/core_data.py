from random import random

import ez_yaml

import pytransit.basics.csv as csv
from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.basics.misc import flatten_once, no_duplicates

class SessionData(LazyDict):
    annotation_path = "/Users/jeffhykin/repos/transit/src/pytransit/genomes/H37Rv_dev.prot_table" # FIXME: default value is for debugging only
    standalone_wigs = []
    combined_wigs = []
    
    @property
    def conditions(self): return no_duplicates(flatten_once(each.conditions for each in self.combined_wigs))
    @property
    def samples(self): return flatten_once(each.samples for each in self.combined_wigs)

universal = LazyDict(
    frame=None,
    session_data=SessionData(),
    busy_running_method=False,
)