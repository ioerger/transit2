# What is this for?

- How to make changes to the pytransit codebase
- What to NOT change
- What tools are available

# Overview

- TLDR; adding a CLI or GUI method
- How does the system work?
    - Where does everything start? Entrypoint
    - What is `global.py`?
    - Whats the difference between specific_tools and generic_tools?
    - Folder Structure
    - How to Manipulate Data (Data structures)
        - CombinedWig
        - Metadata
        - AnnotationFile
        - Gene
        - Genes
- How to make a CLI Method (Tooling)
- Editing a GUI method
- Testing

## TLDR; How to add a new method (CLI or GUI)

1. Copy `src/pytransit/methods/__example__.py` and name it something (I'll call mine `ex1`)
2. Edit `src/pytransit/methods/__init__.py` to import the method, for example I would add `from .ex1 import *` 
3. You *can* (not necessarily should) delete all the example lines and add the following:

```py
from pytransit.globals import logging, cli, gui

# example
@cli.add_command("blah_blah")
def the_blah_blah_command(args, kwargs):
    logging.log(f"you ran the blah_blah command with these args {args}")
```

Running `transit blah_blah` will now print out that message.

NOTE: See the "How to make a CLI Method (Tooling)" section for helper tools, and explaination of args/kwargs.

4. If adding a GUI method
- keep the example code
- search (ctrl/cmd + f) for all the `HANDLE_THIS` comments in the example code
- handle/remove each of the ,HANDLE_THIS comments
- read:
    - the "How to Manipulate Data" section in this guide
    - "Editing a GUI method" section
- done!

## How does the system work?

### Where does everything start? `__main__.py` Explained

All pytransit tasks (gui or cli) start with `__main__.py`
The main things that get run are:
0. `import pytransit.globals` e.g. the globals.py code runs first
1. `import pytransit.transit_cli` then the transit cli code
2. `run_main()` then the run_main function (which is a small wrapper)
3. `clean_args(sys.argv[1:])` how args are handles
4. `main(args, kwargs)` then the actual main function is called

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
```

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

### Generic VS Specific Tools

Generic tools are any python helpers that would be useful on a project unrelated to transit.
For example, removing duplicates from a list (while preserving order) is not specific to transit, and to keep things conceptually simple its very good seperate transit-logic from generic-logic. Which is why something like reading CSV files is inside of `generic_tools/csv.py`.

In contrast, something like parsing a `.wig` is highly specific to transit. Usually specific tools will utilize the generic tools, and then add in logic specific to transit. Thats why something like `CombinedWig.load` is inside of `specific_tools/tnseq_tools.py`

### What is the basic Folder Structure?

The most important part of the folder structure for a hello world will be this:

```sh
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

We've covered `__main__.py`, `globals.py`, `transit_cli.py`, `transit_gui.py` in the sections above, and "specific_tools VS generic_tools" was given its own section above, so `methods/` is the last thing we need to cover for a hello_world.

For a proper `hello_world` command we need to do two things.
1. Create `methods/hello_world.py`
2. Add `from .hello_world import *` to the `methods/__init__.py` file

Once that is done, the command `python ./src/transit.py hello_world jeff` should print out `Hello jeff`

There's a bit more to the folder structure, but we'll cover that after we talk about data and GUI tooling.

## How to Manipulate Data: The Data Structures

After skimming the generic structures, if you have input data (ex: Gff file) structure feel free to skip to particular sections.

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
glycerol_comwig    = combined_wig_object.with_only(condition_names=["Glycerol"])
cholesterol_comwig = combined_wig_object.with_only(condition_names=["Cholesterol"])

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

### AnnotationFile (Gff & Prot Table)

#### How to create the object

```py
annotation = AnnotationFile(path="./somewhere.gff3")
annotation = AnnotationFile(path="./somewhere.prot_table")
```

This class looks at the file extension, and if its neither gff or prot_table then it throws an error explaining what extension is required.

#### Available attributes and methods

```py
# the original filepath 
annotation.path 

# getting descriptions
description = annotation.gene_description(orf_id="Rv0001", fallback_value="-")
# NOTE: if ORF doesnt exist there is NO ERROR, just fallback_value
# if fallback_value is not given the default value is None rather than a string

# 
# a dictionary of named-lists 
# 
annotation.orf_to_info
# example access:
annotation.orf_to_info["Rv0001"].name
annotation.orf_to_info["Rv0001"].description
annotation.orf_to_info["Rv0001"].start_coordinate
annotation.orf_to_info["Rv0001"].end_coordinate
annotation.orf_to_info["Rv0001"].strand
# example 2 access:
name, description, start_coordinate, end_coordinate, strand = annotation.orf_to_info["Rv0001"]

# 
# a list of dictionaries 
# 
# unnecessary? yes, it would be nice to remove this but it
#              was created when consolidating duplicated code
annotation.as_list_of_dicts
# each dictionary has the following keys:
#     "start"
#     "end"
#     "rv"
#     "gene"
#     "strand"
```

### Gene

TODO:

### Genes

TODO:

## How to make a CLI Method (Tooling)

### Global Tools

```py
from pytransit.globals import logging, cli, gui, debugging_enabled, root_folder

logging.log("Howdy")     # use instead of print: it updates GUI if GUI is active
logging.warn("Message")  # does not throw an error, uses stderr in CLI
logging.error("Message") # will exit CLI. In GUI it only throws an error (doesn't crash GUI interface)

gui.is_active     # boolean
root_folder       # absolute path to root folder of transit
debugging_enabled # boolean that is usually used for enabling verbose output
```

Please use `logging.log` instead of `print`

### Argument Tools

What is the difference between `transit_cli.py` and `specific_tools/console_tools.py`?<br>Both `transit_gui.py` and `transit_cli.py` use `console_tools.py`. And `transit_gui.py` and `transit_cli.py` will never be running at the same time (either CLI or GUI)

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

Lets say we want 'hello_world' to be the first command shown whenever the list of of all commands is printed.
The default order is the order they are imported (mostly inside the `methods/__init__.py` file)

However, there is a more direct way to control the order, and to use that we would up the `globals.py` file, and find this part:

```py
@singleton
class cli:
    # the order subcommands are shown can be defined here
    subcommands = {
        ("help",): "This should (and likely is) be replaced by a function elsewhere in the code"
    }
```

Editing it to look something like this will make the `hello_world` command be listed first.

```py
@singleton
class cli:
    # the order subcommands are shown can be defined here
    subcommands = {
        ("hello_world",): "This should (and likely is) be replaced by a function elsewhere in the code",
        ("help",): "This should (and likely is) be replaced by a function elsewhere in the code",
    }
```

If we had subcommands, such as `export hello_world` then we would could change there order like this:

```py
@singleton
class cli:
    # the order subcommands are shown can be defined here
    subcommands = {
        ("hello_world",): "This should (and likely is) be replaced by a function elsewhere in the code",
        ("export", "hello_world",): "This should (and likely is) be replaced by a function elsewhere in the code",
        ("help",): "This should (and likely is) be replaced by a function elsewhere in the code",
    }
```

NOTE: This does mean, if a command is renamed, you'll have to remember to rename it here as well.

## Editing a GUI method

Lets start with trying to edit `methods/resampling.py`

Ideally each analysis method has 4 member functions:
- `output()` this does all the computation, its the main function and usually writes to file
- `from_args()` converts CLI args, then calls the `ouput()` function using those args
- `define_panel()` creates/defines the GUI components
- `from_gui()` grabs data from the GUI, then calls the `ouput()` function using that data

Sometimes there is also a `compute()` which is reserved for pure computation (does not write any output files). This is especially useful for something such as normalization which so that other methods can use the `compute` function without creating side effects.

All of these are not programmatically necessary, but are a helpful convention.

Here's what a skeleton might look like

```py
@misc.singleton
class Method:
    
    # Just like the hello world tutorial 
    @staticmethod
    @cli.add_command("cli_name")
    def from_args(args, kwargs):
        Method.output(
            annotation_path=args[0],
            should_normalize=kwargs["normalize"],
        )
    
    # What happens when menu method is clicked
    @gui.add_menu("Method", "menu_name")
    def define_panel(event):
        ...
        # IMPORTANT: connects "from_gui" to the run button
        panel_helpers.create_run_button(panel, main_sizer, from_gui_function=Method.from_gui)
    
    # What happens when run button is clicked
    @staticmethod
    def from_gui(frame):
        # use panel_helpers to grab data 
        generic_args = GUI_SPECIFIC_STUFF()
        Method.output(generic_args)
        
    # the core function (writes to an output file)
    def output(*, arg1=None, arg2=None):
        # (seems round-about but I promise this style of setting defaults is not frivolous )
        arg1 = arg1 if arg1 is not None else "default value"
        arg2 = arg2 if arg2 is not None else "default value"
        ...
```

TODO: convention:
    from_args
    Run
    from_gui
        panel helpers

## Testing


### Running tests
- run `commands/test/all` to simply test everything
- run `commands/test/all resampling` to test just the CLI tests of resampling
- git add and git commit 

### Adding cli tests
- lets say we're adding one to `resampling`
- copy an example such as `tests/cli_tests/resampling/0001_basic.sh`
- edit it to be the CLI you want to try
- thats

Note: the `|` in `|heapmap` is a hack to get heatmap cli tests to run first