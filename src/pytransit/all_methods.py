from os.path import dirname, basename, isfile
from os import listdir

from super_map import LazyDict

all_methods = LazyDict()
method_names = [basename(name)[:-3] for name in listdir(dirname(__file__)) if isfile(dirname(__file__)+"/"+name) and name not in ('__init__.py')]
for each in method_names:
    exec(f"import pytransit.method.{each} as {each}; all_methods.{each} = {each}")