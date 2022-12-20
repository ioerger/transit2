# What is this for?

This is a document about how to make changes to the pytransit codebase.

Read this before doing tasks such as:

- adding a new command line option
- adding a new input to a method
- creating a new preprocessing method
- adding a new analysis method entirely
- adding a new test for an analysis method

# Overview

- Creating a new CLI method
- Global Interface
- Folder Structure
- Data structures
- Editing a GUI method
- Creating a new GUI method
- Testing

## Adding a new command line option

Adding ONE new command is easy to hack in.<br>But adding many (while staying sane) requires some design, which is what this is going to cover.

Lets start with the **folder structure**. Everything I talk about will be within `src/pytransit/`.

### The Entrypoint: `__main__.py`

All pytransit tasks (gui or cli) start with `__main__.py`
The main things that get run are:
0. `import pytransit.globals`
1. `import pytransit.transit_cli`
2. `run_main()`
3. `clean_args(sys.argv[1:])`
4. `main(args, kwargs)`

The `__main__.py` will end up triggering one of two things:
- A CLI command (just a regular python function)
- Or it will use the `transit_gui.py` file to create a GUI using `TnSeqFrame()`

WHEN SHOULD `__main__.py` BE EDITED:
- (Very rarely in general)
- If adding some alternative to the GUI / CLI interface, such as a web interface
- If the `methods/` folder is renamed
- If the python version check needs to be updated
- If the matplotlib backend needs to be changed

WHEN SHOULD `__main__.py` NOT BE EDITED:
- adding a new CLI command
- changing the order that CLI commands are listed
- updating the default help message


### `globals.py` - Overview

This is the most important file for adding user interfaces to the CLI or GUI.

```py
from pytransit.globals import logging, cli, gui

# example
@cli.add_command("blah_blah")
def the_blah_blah_command(args, kwargs):
    print(f"you ran the blah_blah command with these args {args}")
```

The code above is a fully-functional hello-world for adding a command to transit.
Whenever a user runs `transit blah_blah` it will print out that message.

BUT there are several problems:
- we need to know how `args`, and `kwargs` have been created
- we shouldn't be using `print`
- we need to pick a good file/folder for this code

WHEN SHOULD `globals.py` BE EDITED:
- (Very rarely in general)
- When we want to change the order that commands are listed in `transit help`
- Adding a prefix (like date/time) EVERY time the log function is called
- If there is a new always-avaiable interface (like the status bar in the GUI)
- If global user-settings are added (ex: if we create a "default folder" setting)
- If data needs to be shared between two seperate calls to an analysis methods

WHEN SHOULD `globals.py` NOT BE EDITED:
- adding a new CLI command
- updating the default help message

### `globals.py` - Tooling

Avaiable tools:

```py
from pytransit.globals import logging, cli, gui, debugging_enabled, root_folder

logging.log("Howdy")     # use instead of print: it updates GUI if GUI is active
logging.warn("Message")  # does not throw an error, uses stderr in CLI
logging.error("Message") # will exit CLI. In GUI it only throws an error (doesn't crash GUI interface)

gui.is_active     # boolean
root_folder       # absolute path to root folder of transit
debugging_enabled # boolean that is usually used for enabling verbose output
```

Please use `logging.log` instead of `print`. 

Updated CLI example:

```py
from pytransit.globals import logging, cli, gui

# example
@cli.add_command("blah_blah")
def the_blah_blah_command(args, kwargs):
    logging.log(f"you ran the blah_blah command with these args {args}")
```

### CLI tooling: Arguments

What is the difference between `transit_cli.py` and `specific_tools/console_tools.py`?<br>Conceptually EITHER `transit_cli.py` OR `transit_gui.py` is the star of the show.<br>In contrast `console_tools.py` can be used by both `transit_gui.py` and `transit_cli.py`.<br>In practice, the line is somewhat blurry.

The `args` and `kwargs` values are created by the `clean_args()` function.
All arguments that look like an int are converted to an int.
All arguments that look like an float are converted to an float.
(despite this you may see unnecessary conversions like `int(args[0])` in the codebase)

#### args
- Is a python list
- It will not contain the name of the program (e.g. `arg[0]` is the first normal/regular arugment)
- All keyword args (and their values) will NOT be inside this list

#### kwargs
- `kwargs` is a dictionary, but instead of throwing an error, it will return `None` if given a key that doesn't exist. This makes boolean checks a lot easier.
- There are two types of keyword arguments:
    - `--flag` a flag has two dashes and is always a boolean value
    - `-full-arg value` a "full" keyword argument has only one dash is a key-value pair. It always consumes the following argument. If there is no next argument `clean_args` will throw an error explaining the problem.
- For a keword argument like `--help` there are multiple valid keys. For example: `kwargs["help"]`, `kwargs["-help"]`, `kwargs["--help"]` will all be the same value. NOTE: as consquence, using the same name for both a flag and full keyword (e.g. `--kwarg` and a `-kwarg` in one commmand) does not work.
- EDGECASE: As a consequence of using a dictionary, if a command repeats a full keyword argument, no error will be thrown. The last version of the kwarg will overwrite the other ones.

### Generic VS Specific Tools

Generic tools are any python helpers that would be useful on a project unrelated to transit.
For example, removing duplicates from a list (while preserving order) is not specific to transit, and to keep things conceptually simple its very good seperate transit-logic from generic-logic. Which is why something like reading CSV files is inside of `generic_tools/csv.py`.

In contrast, something like parsing a `.wig` is highly specific to transit. Usually specific tools will utilize the generic tools, and then add in logic specific to transit. Thats why something like `CombinedWig.load` is inside of `specific_tools/tnseq_tools.py`

### CLI tooling: `console_tools.py`

Basically all transit CLI functions should use these helpers:
- `handle_help_flag()`
- `handle_unrecognized_flags()`
- `enforce_number_of_args()`

Lets add them to a hello world example:

```py
from pytransit.globals import logging, cli, gui
from pytransit.specific_tools import console_tools

# what flags are allowed
valid_cli_flags = [ "--excited", ]
# the help message
usage_string = f"""
    {console_tools.subcommand_prefix} hello_world <name>
"""
@cli.add_command("hello_world")
def the_hello_world_command(args, kwargs):
    # checks
    console_tools.handle_help_flag(kwargs, usage_string)
    console_tools.handle_unrecognized_flags(valid_cli_flags, kwargs, usage_string)
    console_tools.enforce_number_of_args(args, usage_string, at_least=1)
    
    # actual hello world
    name = args[0]
    logging.log(f"Hello {name}")
```

### CLI tooling: Sub Commands

Creating sub command is trivial, simply add another string

```py
@cli.add_command("export", "hello_world")
@cli.add_command("hello_world")
def the_hello_world_command(args, kwargs):
    name      = args[0]
    file_path = kwargs["path"]
    # write to file
    with open(file_path, 'w') as the_file:
        the_file.write(f"Hello {name}")
```

This will work for both:
- `python ./src/transit.py hello_world        jeff -path file.txt`
- `python ./src/transit.py export hello_world jeff -path file.txt`

### CLI tooling: The Order of Listed Commands

TODO: example of globals.py

### Running a Hello World: Folder Structure

The most important part of the folder structure will look like this:

```
__main__.py
globals.py
transit_cli.py
transit_gui.py
specific_tools/
generic_tools/
methods/
    __init__.py
    hello_world.py
```

For a proper `hello_world` command we need to do two things.
1. Create `methods/hello_world.py`
2. Add `from .hello_world import *` to the `methods/__init__.py` file

Once that is done the command `python ./src/transit.py hello_world jeff` should print out `Hello jeff`

## How to Manipulate Data: The Data Structures

#### Two generic structures: `LazyDict` and `named_list`

`LazyDict` inherits from `dict` and is fully backwards compatible with all dict methods
`named_list`'s inhert from `list` and are fully backwards compatible with all list methods

`LazyDict`: The only difference from `dict` is `a_dict["thing"] == a_dict.thing`. Its literally only used for simplifying syntax.

`named_list` is slightly different, `named_list()` returns a CLASS.
It is most easily explained with an example:

```py
class Point(named_list(["x", "y"])):
    pass

# example: not all values need to be named
three_d_point = Point([1,2,3])

print(three_d_point[0])  # 1
print(three_d_point[1])  # 2
print(three_d_point[2])  # 3

print(three_d_point.x)   # 1
print(three_d_point.y)   # 2

print(three_d_point)     # [ x=1, y=2, 3 ]
```

### Main Data Structures

Transit has 5 main data structures, all are inside `tnseq_tools.py`.

```
CombinedWig
CombinedWigMetadata
AnnotationFile
Gene
Genes
```

There are some smaller structures like `GffFile` and `Wig`, but they follow from the main data structures

### CombinedWig

Any time you want to work with `.comwig` files, this class is designed to handle all the boilerplate work, such as:
- normalizing
- filtering out conditions
- filtering out wig_ids
- summing across conditions
- getting insertions as a list of rows
- getting insertions as a list of columns
- getting a list of wig_fingerprints
- getting insertions grouped by genes
- getting insertions as a numpy array
- etc

Here's how a combined wig is created

```py
combined_wig_object = CombinedWig(
    main_path="a/file/path.comwig",
    metadata_path="optional/file/path.metadata",
    annotation_path="optional/file/path.tsv",
)
```

Here is reference of the available properties/methods:

```py
combined_wig_object.as_tuple         # (numpy.array(sites), numpy.array(counts_by_wig), wig_fingerprints)
combined_wig_object.rows             # equivalent to the CSV rows of .comwig file; a list of lists, can contain numbers and strings
combined_wig_object.wig_ids          # same order as columns/wig_fingerprints
combined_wig_object.wig_fingerprints # same order as #File: columns
combined_wig_object.read_counts_array[row_index, wig_index]
combined_wig_object.main_path
combined_wig_object.metadata_path   # to get all these it would be [ each.metadata_path for each in gui.combined_wigs ]
combined_wig_object.samples         # list of Wig objects
combined_wig_object.metadata        # CombinedWigMetadata object
combined_wig_object.copy()
combined_wig_object.with_only(condition_names=["name"])
combined_wig_object.with_only(wig_fingerprints=["name"])
combined_wig_object.with_only(wig_ids=["name"])
combined_wig_object.normalized_with(kind="TTR")
combined_wig_object.with_loess_correction()
# TODO: should probably add a .trim(n_terminus=0.0, c_terminus=0.0,) method to CombinedWig
combined_wig_object.averaged(by_conditions=True)
combined_wig_object.averaged(by_genes=True, n_terminus=0, c_terminus=0)
combined_wig_object.summed(by_conditions=True)
combined_wig_object.get_genes( # returns a Genes() object (explained in the Gene's section)
    ignore_codon=True,
    n_terminus=0.0,
    c_terminus=0.0,
    include_nc=False,
    reps="All",
    minread=1,
    genome="",
    transposon="himar1",
)
```

All methods return a copy (e.g. the comwig object is treated as if it is immutable). This makes it very easy to give a copy to a helper function. For example: 

```py
glycerol_comwig    = combined_wig_object.with_only(condition_names=["glycerol"])
cholesterol_comwig = combined_wig_object.with_only(condition_names=["cholesterol"])

from pytransit.specific_tools.transit_tools import calc_gene_means
glycerol_means   , _, _ = calc_gene_means(combined_wig=glycerol_comwig)
cholesterol_means, _, _ = calc_gene_means(combined_wig=cholesterol_comwig)
```

Another example:
```py
ttr_comwig      = combined_wig_object.normalized_with(kind="TTR")
betageom_comwig = combined_wig_object.normalized_with(kind="betageom")

from pytransit.specific_tools.transit_tools import calc_gene_means
ttr_means     , _, _ = calc_gene_means(combined_wig=ttr_comwig)
betageom_means, _, _ = calc_gene_means(combined_wig=betageom_comwig)
```

Another example:
```py
corrected_comwig     = combined_wig_object.with_loess_correction()
not_corrected_comwig = combined_wig_object.copy()

from pytransit.specific_tools.transit_tools import calc_gene_means
corrected_means     , _, _ = calc_gene_means(combined_wig=corrected_comwig)
not_corrected_means, _, _ = calc_gene_means(combined_wig=not_corrected_comwig)
```

Another example:
```py
comwig_of_one_and_two = combined_wig_object.with_only(wig_fingerprints=["path/to/1.comwig", "path/to/2.comwig"])
comwig_of_three       = combined_wig_object.with_only(wig_fingerprints=["path/to/3.comwig", ])

from pytransit.specific_tools.transit_tools import calc_gene_means
one_and_two_means, _, _ = calc_gene_means(combined_wig=comwig_of_one_and_two)
three_means      , _, _ = calc_gene_means(combined_wig=comwig_of_three)
```

CombinedWig objects are desiged to have a Metadata object if the paths is available. Whenever an operation is done (such as `.with_only()` the metadata is also filtered to only include what is mentioned).

### Metadata

Here's how to create a metadata object:

```py
metadata_obj = CombinedWigMetadata(path="path/to/file.metadata")
```

Here is reference of the available properties/methods:

```py
metadata_obj.path             # string or None (is None after .with_only() is called)
metadata_obj.headers          # a list of strings that always (at least) includes:
                              # [ "Condition", "Filename", "Id" ]
metadata_obj.rows             # a list of named lists
metadata_obj.conditions       # a list of Condition objects
metadata_obj.condition_names  # a list of strings
metadata_obj.wig_ids          # a list of strings
metadata_obj.wig_fingerprints # a list of strings
metadata_obj.with_only(condition_names=[])           # returns a CombinedWigMetadata obj
metadata_obj.with_only(wig_fingerprints=[])          # returns a CombinedWigMetadata obj
metadata_obj.condition_names_for(wig_fingerprint="") # returns a list of strings
metadata_obj.condition_names_for(wig_id="")          # returns a list of strings
metadata_obj.id_for(wig_fingerprint="")              # returns a string
metadata_obj.fingerprints_for(condition_name="")     # returns a list of strings
```

### AnnotationFile

TODO:

### Gene

TODO:

### Genes

TODO:

## Editing a GUI method

TODO: convention:
    from_args
    Run
    from_gui
        panel helpers

## Creating a new GUI method

TODO:
- adding a menu item
- result files

## Testing