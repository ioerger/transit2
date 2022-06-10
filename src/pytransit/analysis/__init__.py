# __all__ = []

from os.path import dirname, basename, isfile
from os import listdir
import glob

analysis_names = [basename(name)[:-3] for name in listdir(dirname(__file__)) if isfile(dirname(__file__)+"/"+name) and name not in ('__init__.py', 'base.py')]
__all__        = [basename(path)[:-3] for path in glob.glob(dirname(__file__) + "/*.py") if isfile(path)]

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