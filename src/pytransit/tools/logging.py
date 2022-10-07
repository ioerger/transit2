from pytransit.globals import gui, cli, root_folder, debugging_enabled

def log(message, *args, **kwargs):
    import inspect
    import os
    if gui.is_active:
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
    if gui.is_active:
        set_status(message)

def warn(*args, **kwargs):
    if gui.is_active:
        from pytransit.tools.gui_tools import set_status
    else:
        def set_status(*args, **kwargs): pass
        
    import warnings
    string_output = _print_to_string(*args, **kwargs)
    warnings.warn("\n"+string_output, stacklevel=2)
    last_line = string_output.strip().split("\n")[-1]
    set_status(last_line)

def error(*args, no_traceback=False, **kwargs):
    import traceback
    if gui.is_active:
        from pytransit.tools.gui_tools import set_status
    else:
        def set_status(*args, **kwargs): pass
    
    # this case is different from below because it uses stderr
    if no_traceback and not gui.is_active:
        import sys
        print(*args, file=sys.stderr, **kwargs)
        exit(1)
    else:
        error_message = _print_to_string(*args, **kwargs)
        last_line = error_message.strip().split("\n")[-1]
        try:
            # this looks dumb but 'raise' is needed to get the traceback
            raise TransitError(error_message)
        except Exception as error:
            if not no_traceback:
                traceback.print_exc()
            
            if gui.is_active:
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
