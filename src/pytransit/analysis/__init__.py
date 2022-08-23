from os.path import dirname, basename, isfile
from os import listdir

# these have to be imported selectively... 
# (not all src/pytransit/analysis/*.py files are needed)
dont_include = ('__init__.py', 'base.py')
analysis_names = [basename(name)[:-3] for name in listdir(dirname(__file__)) if isfile(dirname(__file__)+"/"+name) and name not in dont_include ]

# import all the analysis files
for each in analysis_names:
    exec(f"from pytransit.analysis import {each}")

export_methods = {
    "norm": norm.Analysis(),
}
methods = {
    each_name: eval(f"{each_name}.Analysis()")
        for each_name in analysis_names
            if each_name not in export_methods
}
