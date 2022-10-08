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
    from pytransit.basics.lazy_dict import stringify
    indent_string = (" "*by) if isinstance(by, int) else by
    string = string if isinstance(string, str) else stringify(string)
    start = indent_string if not ignore_first else ""
    return start + string.replace("\n", "\n"+indent_string)

def singleton(my_class):                                             
    return my_class()

def human_readable_data(obj):
    import ez_yaml
    import json
    ez_yaml.yaml.version = None
    return ez_yaml.to_string(json.loads(json.dumps(obj)))

def pascal_case_with_spaces(string):
    digits = "1234567890"
    valid_word_contents = "1234567890qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM"
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
            new_string += " "+each_character.upper()
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