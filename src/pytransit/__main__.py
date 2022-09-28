import sys
import traceback

from pytransit.tools import console_tools
from pytransit.universal_data import universal

def main(*args, **kwargs):
    # 
    # Check python version
    # 
    if sys.version_info[0] < 3:
        print("TRANSIT v3.0+ requires python3.6+. To use with python2, please install TRANSIT version < 3.0")
        sys.exit(0)
    
    # 
    # parse args
    # 
    if True:
        # 
        # extract debug value
        # 
        DEBUG = "--debug" in sys.argv
        if DEBUG:
            sys.argv.remove("--debug")
            kwargs.pop("-debug")
        
        # 
        # version command
        # 
        if not args and ("v" in kwargs or "-version" in kwargs):
            import pytransit
            print("Version: {0}".format(pytransit.__version__))
            sys.exit(0)
        
        # 
        # help command
        # 
        if not args and ("h" in kwargs or "-help" in kwargs):
            print_help()
            sys.exit(0)

    # 
    # GUI Mode
    # 
    if not (args or kwargs):
        universal.interface = "gui"
        # Tried GUI but no wxPython
        if not console_tools.check_if_has_wx():
            print("Please install wxPython to run in GUI Mode. (pip install wxPython)")
            print("")
            print_help()
        # WX is available
        else:
            import matplotlib
            import matplotlib.pyplot
            import pytransit.transit_gui as transit_gui
            import wx

            matplotlib.use("WXAgg")

            print("Running in GUI Mode")
            app = wx.App(False)

            # create an object of CalcFrame
            frame = transit_gui.TnSeqFrame(None, DEBUG)
            # show the frame
            frame.Show(True)
            frame.Maximize(True)

            # start the applications
            app.MainLoop()
    # 
    # Console mode
    # 
    else:
        universal.interface = "console"
        import matplotlib
        from pytransit.methods.analysis import methods as analysis_methods
        from pytransit.methods.export   import methods as export_methods
        from pytransit.methods.convert  import methods as convert_methods
        
        def run(method, args):
            if hasattr(method, "method"): # TODO: remove this once the convert/export methods have been updated (probably in a few weeks - Oct 13st) --Jeff
                method = method.method
                
            setup_object = None
            try:
                setup_object = method.from_args(args, kwargs)
            # dont show traceback for invalid argument errors
            except console_tools.InvalidArgumentException as error:
                print("Error: %s" % str(error))
                print(method.usage_string)
                sys.exit(1)
            # show traceback and usage for all other kinds of errors
            except Exception as error:
                print("Error: %s" % str(error))
                traceback.print_exc()
                print(method.usage_string)
                sys.exit(1)
            
            if setup_object:
                setup_object.Run()
                sys.exit(0)
        
        def check_if_missing(kind, selected_name, methods):
            if selected_name not in methods:
                print(f"Error: Need to specify the {kind} method.")
                print(f"Please use one of the known methods (or see documentation to add a new one):")
                for each_export_method in methods:
                    print(f"\t - {each_export_method}")
                print(f"Usage: python {sys.argv[0]} {kind} <method>")
                sys.exit(1)
        
        matplotlib.use("Agg")
        method_name, *args = args
        # 
        # Analysis
        # 
        if method_name in analysis_methods:
            run(
                method=analysis_methods[method_name],
                args=args,
            )
        
        # 
        # Export
        # 
        if method_name.lower() == "export":
            export_method_name = ""
            if len(args) >= 1:
                export_method_name, *args = args
            check_if_missing(kind="export", selected_name=export_method_name, methods=export_methods)
            run(
                method=export_methods[export_method_name],
                args=args[1:],  # skip the first argument for some reason
            )
        # 
        # Convert
        # 
        elif method_name.lower() == "convert":
            convert_method_name = ""
            if len(args) >= 1:
                convert_method_name = args[1]
            
            check_if_missing(kind="convert", selected_name=convert_method_name, methods=convert_methods)
            run(
                method=export_methods[export_method_name],
                args=args,
            )
        # 
        # Error
        # 
        else:
            print(f"Error: The '{method_name}' method is unknown.")
            print("Please use one of the known methods (or see documentation to add a new one):")
            for each_method_name in analysis_methods:
                print("\t - %s" % each_method_name)
            print(f"Usage: python {sys.argv[0]} <method>")

def run_main():
    from pytransit.tools.console_tools import clean_args
    (args, kwargs) = clean_args(sys.argv[1:])
    main(*args, **kwargs)

def print_help():
    from pytransit.methods.analysis import methods as analysis_methods
    print("To run in Console Mode, try 'transit <method>' with one of the following methods:")
    print("Analysis methods: ")
    for each_analysis_method in analysis_methods:
        ## TODO :: Move normalize to separate subcommand?
        if each_analysis_method == "normalize":
            continue
        print(f"    - {each_analysis_method}")
    print("Other methods: ")
    print("    - normalize")
    print("    - convert")
    print("    - export")
    print(f"Usage: python {sys.argv[0]} <method>")
    print(f"Example: python {sys.argv[0]} normalize --help")
        
if __name__ == "__main__":
    main()
