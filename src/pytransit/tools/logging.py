from pytransit.universal_data import universal

def log(message, *args, **kwargs):
    import inspect
    import os
    if universal.interface == 'gui':
        from pytransit.tools.gui_tools import set_status
    else:
        def set_status(*args, **kwargs): pass
        
    
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
        set_status(message)

def warn(*args, **kwargs):
    if universal.interface == 'gui':
        from pytransit.tools.gui_tools import set_status
    else:
        def set_status(*args, **kwargs): pass
        
    import warnings
    string_output = _print_to_string(*args, **kwargs)
    warnings.warn(string_output)
    last_line = string_output.strip().split("\n")[-1]
    set_status(last_line)

def error(*args, **kwargs):
    import traceback
    if universal.interface == 'gui':
        from pytransit.tools.gui_tools import set_status
    else:
        def set_status(*args, **kwargs): pass
    
    error_message = _print_to_string(*args, **kwargs)
    last_line = string_output.strip().split("\n")[-1]
    try:
        # this looks dumb but 'raise' is needed to get the traceback
        raise TransitError(error_message)
    except Exception as error:
        traceback.print_exc()
        if universal.interface == 'gui':
            set_status(last_line)
        else:
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
