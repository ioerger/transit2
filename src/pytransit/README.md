# What is this for?

This is a document about how to make changes to the pytransit codebase.

Read this before doing any of the following:

- adding a new command line option
- adding a new input to a method
- creating a new preprocessing method
- adding a new analysis method entirely

## Adding a new command line option

Adding ONE new option is easy. <br> Adding many options (while staying sane) requires some forethought.

Lets start with the **folder structure**. Everything I talk about will be within `src/pytransit/`.

### The Entrypoint `__main__.py`

All pytransit tasks start with `__main__.py`
The function/file sequence looks like:
0. `import pytransit.globals`
1. `import pytransit.transit_cli`
2. `run_main()`
3. `clean_args(sys.argv[1:])`
4. `main(args, kwargs)`

The `__main__.py` will trigger one of two things:
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


### The `globals.py`

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
- If there is a new always-avaiable interface (like the status bar in the GUI)
- If global user-settings are added (ex: if we create a "default folder" setting)
- If data needs to be shared between two seperate calls to an analysis methods

WHEN SHOULD `globals.py` NOT BE EDITED:
- adding a new CLI command
- updating the default help message


### CLI tool


- Folder structure
    - `globals.py`
    - `data/`
    - `specific_tools/`
    - `gen`
    - `methods/`
- Tools
    - Logging
    - Testing
    - General Tools
        - LazyDict
        - named_list
        - csv.read
        - csv.write
        - no_duplicates
    - Specific Tools
        - console_tools.enforce_number_of_args
        - console_tools.handle_unrecognized_flags
        - 
- Control Flow
- Data Structures
- Analysis Method Structure


General explanation
Folder structure
General tools
Specific tools
Method template
Logging
TODO: move progress_update inside of logging
Data access
Global data
Debugging flag
root_folder
Cli object
add_command
Gui object
add_menu
CombinedWig
filtering
normalizing 
summing
Metadata
Annotation
Genes
Gene
Method Template
Inputs as a dictionary
file_kind=Method.identifier
from_args
console_tools.handle_help_flag, etc
from_gui
panel_helpers
Run
Result file format / usage
