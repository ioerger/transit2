def indent(string, by, ignore_first=False):
    indent_string = (" "*by) if isinstance(by, int) else by
    string = string if isinstance(string, str) else stringify(string)
    start = indent_string if not ignore_first else ""
    return start + string.replace("\n", "\n"+indent_string)

def stringify(value):
    class Map:
        pass
    
    onelineify_threshold = 50 # characters (of inner content)
    length = 0
    if isinstance(value, str):
        return f'"{value}"'
    elif isinstance(value, Map):
        if len(value) == 0:
            return "{}"
        items = value if isinstance(value, Map) else value.items()
        output = "{\n"
        for each_key, each_value in items:
            element_string = stringify(each_key) + ": " + stringify(each_value)
            length += len(element_string)+2
            output += indent(element_string, by=4) + ", \n"
        output += "}"
        if length < onelineify_threshold:
            output = output.replace("\n    ","").replace("\n","")
        return output
    elif isinstance(value, dict):
        if len(value) == 0:
            return "{}"
        items = value if isinstance(value, Map) else value.items()
        output = "{\n"
        for each_key, each_value in items:
            element_string = stringify(each_key) + ": " + stringify(each_value)
            length += len(element_string)+2
            output += indent(element_string, by=4) + ", \n"
        output += "}"
        if length < onelineify_threshold:
            output = output.replace("\n    ","").replace("\n","")
        return output
    elif isinstance(value, list):
        if len(value) == 0:
            return "[]"
        output = "[\n"
        for each_value in value:
            element_string = stringify(each_value)
            length += len(element_string)+2
            output += indent(element_string, by=4) + ", \n"
        output += "]"
        if length < onelineify_threshold:
            output = output.replace("\n    ","").replace("\n","")
        return output
    elif isinstance(value, set):
        if len(value) == 0:
            return "set([])"
        output = "set([\n"
        for each_value in value:
            element_string = stringify(each_value)
            length += len(element_string)+2
            output += indent(element_string, by=4) + ", \n"
        output += "])"
        if length < onelineify_threshold:
            output = output.replace("\n    ","").replace("\n","")
        return output
    elif isinstance(value, tuple):
        if len(value) == 0:
            return "tuple()"
        output = "(\n"
        for each_value in value:
            element_string = stringify(each_value)
            length += len(element_string)+2
            output += indent(element_string, by=4) + ", \n"
        output += ")"
        if length < onelineify_threshold:
            output = output.replace("\n    ","").replace("\n","")
        return output
    else:
        try:
            debug_string = value.__repr__()
        except Exception as error:
            from io import StringIO
            import builtins
            string_stream = StringIO()
            builtins.print(value, file=string_stream)
            debug_string = string_stream.getvalue()
        
        # TODO: handle "<slot wrapper '__repr__' of 'object' objects>"
        if debug_string.startswith("<class") and hasattr(value, "__name__"):
            return value.__name__
        if debug_string.startswith("<function <lambda>"):
            return "(lambda)"
        if debug_string.startswith("<function") and hasattr(value, "__name__"):
            return value.__name__
        if debug_string.startswith("<module") and hasattr(value, "__name__"):
            _, *file_path, _, _ = debug_string.split(" ")[-1]
            file_path = "".join(file_path)
            return f"module(name='{value.__name__}', path='{file_path}')"
        
        space_split = debug_string.split(" ")
        if len(space_split) >= 4 and debug_string[0] == "<" and debug_string[-1] == ">":
            
            if space_split[-1].startswith("0x") and space_split[-1] == "at":
                _, *name_pieces = space_split[0]
                *parts, name = "".join(name_pieces).split(".")
                parts_str = ".".join(parts)
                return f'{name}(from="{parts_str}")'
        
        return debug_string

defaulters = {}
class LazyDict(dict):
    
    def __init__(self, *args, **kwargs):
        # default vaue
        super(LazyDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
        defaulters[id(self)] = lambda key: None
    
    def __getitem__(self, key):
        # defaulter value
        defaulter = defaulters.get(id(self))
        if defaulter:
            if key not in self.__dict__:
                return defaulter(key)
        return self.__dict__.get(key)
        
    def __delitem__(self, key):
        try:
            del self.__dict__[key]
        except Exception as error:
            pass
    
    def __str__(self):
        if len(self.__dict__) == 0:
            return "{}"
        return stringify(self.__dict__)
    
    def __repr__(self):
        return self.__str__()
    
    def merge(self, other_dict=None, **kwargs):
        other_dict = other_dict or {}
        self.__dict__.update(other_dict)
        self.__dict__.update(kwargs)
        return self
    
    def update(self, other_dict):
        for each_key, each_value in other_dict.items():
            self[each_key] = each_value
        return self
    
    def setdefault(self, *args, **kwargs):
        if len(args) == 1:
            if callable(args[0]):
                defaulters[id(self)] = args[0]
            else:
                defaulters[id(self)] = lambda key: args[0]
            return self
        else:
            return self.__dict__.setdefault(*args, **kwargs)
