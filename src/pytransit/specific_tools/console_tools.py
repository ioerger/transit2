import os
import sys

from pytransit.generic_tools import command_line
probably_used_transit_executable_directly = os.path.basename(sys.argv[0]) == "transit"
if probably_used_transit_executable_directly:
    full_commandline_command = "transit " + " ".join([ f"'{each}'" for each in sys.argv[1:]])
    subcommand_prefix = command_line.color(
        f" transit ",
        foreground="bright_blue"
    )
else:
    full_commandline_command = f"python3 " + " ".join([ f"'{each}'" for each in sys.argv])
    subcommand_prefix = command_line.color(
        f" python3 {sys.argv[0]} ",
        foreground="bright_blue"
    )




class InvalidArgumentException(Exception):
    def __init__(self, *args, **kwargs):
        super(InvalidArgumentException, self).__init__(*args, **kwargs)

def clean_args(rawargs):
    """Returns a list and a dictionary with positional and keyword arguments.

    -This function assumes flags must start with a "-" and and cannot be a 
        number (but can include them).
    
    -Flags should either be followed by the value they want to be associated 
        with (i.e. -p 5) or will be assigned a value of True in the dictionary.

    -The dictionary will map flags to the name given minus ONE "-" sign in
        front. If you use TWO minus signs in the flag name (i.e. --verbose), 
        the dictionary key will be the name with ONE minus sign in front 
        (i.e. {"-verbose":True}).
    

    Arguments:
        rawargs (list): List of positional/keyword arguments. As obtained from
                         sys.argv.

    Returns:
        list: List of positional arguments (i.e. arguments without flags),
                in order provided.
        dict: Dictionary mapping flag (key is flag minus the first "-") and
                their values.

    """
    args = []
    kwargs = {}
    count = 0
    # Loop through list of arguments
    while count < len(rawargs):
        # If the current argument starts with "-", then it's probably a flag
        if rawargs[count].startswith("-"):
            # Check if next argument is a number
            try:
                temp = float(rawargs[count + 1])
                next_is_number = True
            except:
                next_is_number = False

            still_not_finished = count + 1 < len(rawargs)
            if still_not_finished:
                next_is_not_argument = not rawargs[count + 1].startswith("-")
                next_looks_like_list = len(rawargs[count + 1].split(" ")) > 1
            else:
                next_is_not_argument = True
                next_looks_like_list = False

            # If still things in list, and they look like arguments to a flag, add them to dict
            kwarg_name = rawargs[count][1:]
            if still_not_finished and (
                next_is_not_argument or next_looks_like_list or next_is_number
            ):
                kwargs[kwarg_name] = rawargs[count + 1]
                count += 1
            # Else it's a flag but without arguments/values so assign it True
            else:
                kwargs[kwarg_name] = True
        # Else, it's probably a positional arguement without flags
        else:
            args.append(rawargs[count])
        count += 1
    
    # make the --name vs -name irrelevent 
    for each_key, each_value in list(kwargs.items()):
        if each_key.startswith("--"):
            kwargs[each_key[1:]] = each_value
            kwargs[each_key[2:]] = each_value
        elif each_key.startswith("-"):
            kwargs[each_key[1:]] = each_value
            kwargs["-"+each_key] = each_value
        else:
            kwargs["--"+each_key] = each_value
            kwargs["-"+each_key] = each_value
            
    return (args, kwargs)

def handle_help_flag(kwargs, usage_string):
    if kwargs.get("-help", False):
        print(usage_string)
        exit(0)

def enforce_number_of_args(args, usage_string, *, exactly=None, at_least=type(None)):
    number_of_args = len(args)
    if at_least != type(None):
        if number_of_args < at_least:
            print(f"Wrong number of arguments ({number_of_args} < {at_least}), see usage below\n")
            print(usage_string)
            exit(1)
    else:
        if number_of_args != exactly:
            print(f"Wrong number of arguments ({number_of_args} != {exactly}), see usage below\n")
            print(usage_string)
            exit(1)

def handle_unrecognized_flags(flags, kwargs, usage_string):
    possible_flags = list(flags)
    for each_name in flags:
        if each_name.startswith("--"):
            possible_flags.append(each_name[1:])
            possible_flags.append(each_name[2:])
        elif each_name.startswith("-"):
            possible_flags.append(each_name[1:])
            possible_flags.append("-"+each_name)
        else:
            possible_flags.append("--"+each_name)
            possible_flags.append("-"+each_name)
    
    for flag_string in kwargs.keys():
        if flag_string not in possible_flags:
            closest_match, *_            = misc.levenshtein_distance_sort(word=flag_string, other_words=possible_flags)
            print(f"{flag_string} was not one of the available flags: {' '.join(flags)}")
            raise Exception(f'''unrecognized flag: {flag_string}\nMaybe you meant: {closest_match}\n\n{usage_string}''')

def check_if_has_wx():
    """
    This is designed to be lazy and more lightweight than the HAS_WX variable inside of transit_tools
    """
    cache = getattr(check_if_has_wx, "cache", None)
    if type(cache) == type(None):
        try:
            import wx
            version = int(wx.version()[0]) # attempt to get version
            check_if_has_wx.cache = True
        except Exception as error:
            check_if_has_wx.cache = False
    cache = check_if_has_wx.cache
    return cache
