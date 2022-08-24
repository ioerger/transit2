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

def indent(string, by, ignore_first=False):
    from pytransit.basics.lazy_dict import stringify
    indent_string = (" "*by) if isinstance(by, int) else by
    string = string if isinstance(string, str) else stringify(string)
    start = indent_string if not ignore_first else ""
    return start + string.replace("\n", "\n"+indent_string)
