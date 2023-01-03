#!/usr/bin/env python3
import os
import sys
import traceback

from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit import transit_cli
import pytransit.methods

# this wrapper exist to help with test cases. Might be good to change the test cases to elimate the need for it
def main(args, kwargs):
    from collections import defaultdict
    kwargs = defaultdict(lambda :None, kwargs)
    
    # 
    # Check python version
    # 
    if sys.version_info[0] < 3:
        print("TRANSIT v3.0+ requires python3.6+. To use with python2, please install TRANSIT version < 3.0")
        sys.exit(0)
    
    # 
    # parse global args
    # 
    transit_cli.parse_global_flags(args, kwargs)

    # 
    # GUI Mode
    # 
    if not args:
        import pytransit.specific_tools.console_tools as console_tools
        # Tried GUI but no wxPython
        if not console_tools.check_if_has_wx():
            print("Please install wxPython to run in GUI Mode. (pip install wxPython)")
            print("")
            print("To run in Console Mode, try 'transit <method>' with one of the following methods:")
            transit_cli.help_command(args, kwargs)
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
            frame = transit_gui.TnSeqFrame(None)
            # show the frame
            frame.Show(True)
            frame.Maximize(True)

            # start the applications
            app.MainLoop()
    # 
    # Console mode
    # 
    else:
        import matplotlib
        matplotlib.use("Agg")
        
        # try to match a subcommand
        transit_cli.find_and_run_subcommand(args, kwargs)
        
        # runtime only gets here if no subcommands were matched
        transit_cli.help_command(args, kwargs)

def run_main():
    from pytransit.specific_tools.console_tools import clean_args
    (args, kwargs) = clean_args(sys.argv[1:])
    main(args, kwargs)

if __name__ == "__main__":
    run_main()