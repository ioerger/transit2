def open_url(url):
    import sys
    import os
    import subprocess
    
    output = ""
    error = ""
    try:
        if sys.platform.startswith("darwin"): # OSX
            args = ["open", url]
            output, error = subprocess.Popen(
                args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
        elif os.name == "nt": # Windows
            os.startfile(url)
        elif os.name == "posix": # Linux
            args = ["xdg-open", url]
            output, error = subprocess.Popen(
                args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if "not found" in error:
                args = ["exo-open", url]
                output, error = subprocess.Popen(
                    args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()

    except Exception as error:
        raise Exception(f'''Error opening {url} \n{error}''')


def is_iterable(thing):
    # https://stackoverflow.com/questions/1952464/in-python-how-do-i-determine-if-an-object-is-iterable
    try:
        iter(thing)
    except TypeError:
        return False
    else:
        return True

def is_generator_like(thing):
    return is_iterable(thing) and not isinstance(thing, (str, bytes))

def iteratively_flatten_once(items):
    for each in items:
        if is_generator_like(each):
            yield from each
        else:
            yield each

def flatten_once(items):
    output = list(iteratively_flatten_once(items))
    return output

def no_duplicates(items): # preserving order
    copy = []
    for each in items:
        if each in copy:
            continue
        copy.append(each)
    return copy

def line_count_of(file_path):
    # from stack overflow "how to get a line count of a large file cheaply"
    def _make_gen(reader):
        while 1:
            b = reader(2**16)
            if not b: break
            yield b
    with open(file_path, "rb") as file:
        count = sum(buf.count(b"\n") for buf in _make_gen(file.raw.read))
    
    return count

def indent(string, by="    ", ignore_first=False):
    from pytransit.generic_tools.lazy_dict import stringify
    indent_string = (" "*by) if isinstance(by, int) else by
    string = string if isinstance(string, str) else stringify(string)
    start = indent_string if not ignore_first else ""
    return start + string.replace("\n", "\n"+indent_string)

def to_pure(an_object, recursion_help=None):
    # 
    # infinte recursion prevention
    # 
    top_level = False
    if recursion_help is None:
        top_level = True
        recursion_help = {}
    class PlaceHolder:
        def __init__(self, id):
            self.id = id
        def eval(self):
            return recursion_help[key]
    object_id = id(an_object)
    # if we've see this object before
    if object_id in recursion_help:
        # if this value is a placeholder, then it means we found a child that is equal to a parent (or equal to other ancestor/grandparent)
        if isinstance(recursion_help[object_id], PlaceHolder):
            return recursion_help[object_id]
        else:
            # if its not a placeholder, then we already have cached the output
            return recursion_help[object_id]
    # if we havent seen the object before, give it a placeholder while it is being computed
    else:
        recursion_help[object_id] = PlaceHolder(object_id)
    
    parents_of_placeholders = set()
    
    # 
    # optional torch tensor converter
    # 
    if hasattr(an_object, "__class__") and hasattr(an_object.__class__, "__name__"):
        if an_object.__class__.__name__ == "Tensor":
            try:
                import torch
                if isinstance(an_object, torch.Tensor):
                    an_object = an_object.detach().cpu()
            except Exception as error:
                pass
    # 
    # main compute
    # 
    return_value = None
    # base case 1 (iterable but treated like a primitive)
    if isinstance(an_object, str):
        return_value = an_object
    # base case 2 (exists because of scalar numpy/pytorch/tensorflow objects)
    elif hasattr(an_object, "tolist"):
        return_value = an_object.tolist()
    else:
        # base case 3
        if not is_iterable(an_object):
            return_value = an_object
        else:
            if isinstance(an_object, dict):
                return_value = {
                    to_pure(each_key, recursion_help) : to_pure(each_value, recursion_help)
                        for each_key, each_value in an_object.items()
                }
            else:
                return_value = [ to_pure(each, recursion_help) for each in an_object ]
    
    # convert iterables to tuples so they are hashable
    if is_iterable(return_value) and not isinstance(return_value, dict) and not isinstance(return_value, str):
        return_value = tuple(return_value)
    
    # update the cache/log with the real value
    recursion_help[object_id] = return_value
    #
    # handle placeholders
    #
    if is_iterable(return_value):
        # check if this value has any placeholder children
        children = return_value if not isinstance(return_value, dict) else [ *return_value.keys(), *return_value.values() ]
        for each in children:
            if isinstance(each, PlaceHolder):
                parents_of_placeholders.add(return_value)
                break
        # convert all the placeholders into their final values
        if top_level == True:
            for each_parent in parents_of_placeholders:
                iterator = enumerate(each_parent) if not isinstance(each_parent, dict) else each_parent.items()
                for each_key, each_value in iterator:
                    if isinstance(each_parent[each_key], PlaceHolder):
                        each_parent[each_key] = each_parent[each_key].eval()
                    # if the key is a placeholder
                    if isinstance(each_key, PlaceHolder):
                        value = each_parent[each_key]
                        del each_parent[each_key]
                        each_parent[each_key.eval()] = value
    
    # finally return the value
    return return_value

def singleton(my_class):                                             
    return my_class()

def human_readable_data(obj):
    import ez_yaml
    import json
    ez_yaml.yaml.version = None
    try:
        return ez_yaml.to_string(obj)
    except Exception as error:
        return ez_yaml.to_string(to_pure(obj))

def pascal_case_with_spaces(string):
    digits = "1234567890-"
    valid_word_contents = "1234567890qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM-"
    new_string = " "
    # get pairwise elements
    for each_character in string:
        prev_character = new_string[-1]
        prev_is_lowercase = prev_character.lower() == prev_character
        each_is_uppercase = each_character.lower() != each_character
        
        # remove misc characters (handles snake case, kebab case, etc)
        if each_character not in valid_word_contents:
            new_string += " "
        # start of word
        elif prev_character not in valid_word_contents:
            new_string += each_character.upper()
        # start of number
        elif prev_character not in digits and each_character in digits:
            new_string += each_character
        # end of number
        elif prev_character in digits and each_character not in digits:
            new_string += each_character.upper()
        # camel case
        elif prev_is_lowercase and each_is_uppercase:
            new_string += " "+each_character.upper()
        else:
            new_string += each_character
    
    # flatten out all the whitespace
    new_string = new_string.strip()
    while "  " in new_string:
        new_string = new_string.replace("  "," ")
    
    return new_string

def levenshtein_distance_sort(*, word, other_words):
    word = word.lower()
    def levenshtein_distance(s1, s2):
        # https://stackoverflow.com/questions/2460177/edit-distance-in-python
        if len(s1) > len(s2):
            s1, s2 = s2, s1
        
        distances = range(len(s1) + 1)
        for i2, c2 in enumerate(s2):
            distances_ = [i2+1]
            for i1, c1 in enumerate(s1):
                if c1 == c2:
                    distances_.append(distances[i1])
                else:
                    distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
            distances = distances_
        return distances[-1]
    
    prioritized = sorted(other_words, key=lambda each_other: levenshtein_distance(word, each_other))
    return prioritized

def all_equal(a_list):
    if len(a_list) == 0:
        return True
    
    prev = a_list[0]
    for each in a_list[1:]:
        if prev != each:
            return False
        prev = each
    
    return True

def all_different(a_list):
    if len(a_list) == 0:
        return True
    
    prev = a_list[0]
    for each in a_list[1:]:
        if prev == each:
            return False
        prev = each
    
    return True

def remove_common_prefix(list_of_strings):
    def all_equal(a_list):
        if len(a_list) == 0:
            return True
        
        prev = a_list[0]
        for each in a_list:
            if prev != each:
                return False
            prev = each
        
        return True
    
    shortest_path_length = min([ len(each_path) for each_path in list_of_strings ])
    longest_common_path_length = shortest_path_length
    while longest_common_path_length > 0:
        # binary search would be more efficient but its fine
        longest_common_path_length -= 1
        if all_equal([ each[0:longest_common_path_length] for each_path in list_of_strings ]):
            break
    
    return [ each[longest_common_path_length:] for each_path in list_of_strings ]

def merge_dicts(*args, **kwargs):
    new_dict = {}
    for each in args:
        new_dict.update(each)
    new_dict.update(kwargs)
    return new_dict

def str_is_int(string):
    return string.lstrip('+-').isdigit()

def str_is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def invert_dict(existing_dict):
    new_dict = {}
    for each_key, each_value in existing_dict.items():
        if each_value not in new_dict:
            new_dict[each_value] = []
        new_dict[each_value].append(each_key)
    return new_dict

def inject_path_extension(path, *, extension):
    import os
    parent_folders = os.path.dirname(path)
    full_name = os.path.basename(path)
    parts = full_name.split(".")
    parts.insert(1,extension)
    if parent_folders:
        return parent_folders+"/"+".".join(parts)
    else:
        return ".".join(parts)

from pytransit.generic_tools.cool_cache import cache, settings
settings.default_folder = None # use in-memory storage



go = None
try:
    import plotly.graph_objs as go
except Exception as error:
    pass
