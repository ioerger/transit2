from pytransit.universal_data import universal

if universal.interface == 'gui':
    from pytransit.tools.gui_tools import set_status
else:
    def set_status(*args, **kwargs): pass

def log(message, *args, **kwargs):
    import inspect
    import os
    message = f"{message} "+ " ".join([ f"{each}" for each in args])
    
    # get some context as to who is creating the message
    stack             = inspect.stack()
    caller_frame_info = stack[1]
    file_name         = ""
    caller_name       = ""
    try: file_name = os.path.basename(caller_frame_info.filename)
    except Exception as error: pass # sometimes the caller doesn't have a file name (ex: REPL)
    try: caller_name = caller_frame_info.function
    except Exception as error: pass # sometimes the caller doesn't have a function name (ex: lambda)
    
    # remove the .py extension
    if file_name[len(file_name)-3:len(file_name)] == ".py":
        file_name = file_name[0:len(file_name)-3]
    
    print(f'[{file_name}:{caller_name}()]', message, flush=True, **kwargs)
    if universal.interface == 'gui':
        import pytransit.tools.gui_tools as gui_tools
        gui_tools.set_status(message)

def warn(*args, **kwargs):
    import warnings
    string_output = _print_to_string(*args, **kwargs)
    warnings.warn(string_output)
    last_line = string_output.strip().split("\n")[-1]
    set_status(last_line)

def error(*args, **kwargs):
    import traceback
    if universal.interface != 'gui':
        try:
            # this looks dumb but its needed to get the traceback
            raise TransitError(*args, **kwargs)
        except Exception as error:
            traceback.print_exc()
            print(method.usage_string)
            exit(1)

class TransitError(Exception):
    pass
    
def _print_to_string(*args, **kwargs):
    from io import StringIO
    string_stream = StringIO()
    # dump to string
    print(*args, **{ "flush": True, **kwargs, "file":string_stream })
    output_str = string_stream.getvalue()
    string_stream.close()
    return output_str
