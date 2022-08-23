from os.path import dirname, basename, isfile
from os import listdir

# these have to be imported selectively... 
# (not all src/pytransit/analysis/*.py files are needed)
dont_include = ('__init__.py', 'base.py')
analysis_names = [basename(name)[:-3] for name in listdir(dirname(__file__)) if isfile(dirname(__file__)+"/"+name) and name not in dont_include ]

source = "pytransit.analysis"
attribute = ".Analysis()"
class MethodsCollection(dict):
    def __getitem__(self, key):
        exec(f"""
            \nglobal MethodsCollection
            \nfrom {source} import {key}
            \nMethodsCollection.__temp_output__ = {key}{attribute}
        """, globals(), locals())
        super(MethodsCollection, self).__setitem__(key, MethodsCollection.__temp_output__)
        return MethodsCollection.__temp_output__

methods = MethodsCollection({ name : None for name in analysis_names })