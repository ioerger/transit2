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

class RawKwargs: pass
def clean_args(raw_args):
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
        raw_args (list): List of positional/keyword arguments. As obtained from
                         sys.argv.

    Returns:
        list: List of positional arguments (i.e. arguments without flags),
                in order provided.
        dict: Dictionary mapping flag (key is flag minus the first "-") and
                their values.

    """
    from pytransit.generic_tools import misc
    from collections import defaultdict
    
    args = []
    kwargs = defaultdict(lambda *args: None)
    kwargs[RawKwargs] = {}
    count = 0
    
    raw_args = list(raw_args)
    while len(raw_args) > 0:
        next_raw_argument = raw_args.pop(0)
        # If the current argument starts with "-", then it's a key
        if next_raw_argument.startswith('--'):
            if len(raw_args) == 0:
                raise Exception(f'''
                    
                    This argument: {next_raw_argument}
                    expects a value after it (ex: --key value), however it was the last argument
                    Maybe you meant: {next_raw_argument[1:]}
                    (which is just a flag, e.g. no value)
                    
                ''')
            
            next_arg_is_possibly_negative_number = False
            try:
                float(raw_args[0])
                next_arg_is_possibly_negative_number = True
            except Exception as error: pass
            
            if raw_args[0].startswith('-') and not next_arg_is_possibly_negative_number:
                raise Exception(f'''
                    
                    This argument: {next_raw_argument}
                    expects a value after it (-key value)
                    However it was follow by another key/flag (--key {raw_args[0]})
                    Maybe this would be valid: -{next_raw_argument} {raw_args[0]}
                    (which is just a flag, e.g. no value)
                    
                ''')
            # consume the next element as the value
            kwargs[next_raw_argument] = raw_args.pop(0)
        # its a flag
        elif next_raw_argument.startswith("-") and not misc.str_is_float(next_raw_argument):
            kwargs[next_raw_argument] = True
        # Else, it's a positional argument without flags
        else:
            args.append(next_raw_argument)
    
    # 
    # convert number arguments to be numbers
    # 
    for each_index, each_value in enumerate(list(args)):
        if misc.str_is_int(each_value):
            args[each_index] = int(each_value)
        if misc.str_is_float(each_value):
            args[each_index] = float(each_value)
    kwargs[RawKwargs] = dict({**kwargs})
    for each_key, each_value in kwargs.items():
        if isinstance(each_value, str):
            if misc.str_is_int(each_value):
                kwargs[each_key] = int(each_value)
            if misc.str_is_float(each_value):
                kwargs[each_key] = float(each_value)
    
    # 
    # make the -name vs -name irrelevent 
    # 
    for each_key, each_value in list(kwargs.items()):
        if isinstance(each_key, str):
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
    if any([ not each.startswith('-') for each in flags]):
        print(f'''usage_string = {usage_string}''')
        raise Exception(f'''Developer issue: {flags} contains elements that are missing the - or -- prefix''')
    kwargs = kwargs[RawKwargs]
    possible_flags = flags
    for flag_string in kwargs.keys():
        if isinstance(flag_string, str):
            if flag_string not in possible_flags:
                from pytransit.generic_tools import misc
                closest_match, *_ = misc.levenshtein_distance_sort(word=flag_string, other_words=possible_flags)
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

def string_arg_to_list(argument):
    if isinstance(argument, str):
        return argument.split(",")
    else:
        return None