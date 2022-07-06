window = None

def bind_to(wxPythonObj, event):
    def wrapper2(function_to_attach):
        wxPythonObj.Bind(event, function_to_attach)
        return function_to_attach
    return wrapper2


import traceback
def transit_handle_error(error_obj):
    traceback.print_exc()
    # TODO: send error object 
    # error_obj.args