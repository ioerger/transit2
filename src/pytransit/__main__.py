#!/usr/bin/env python3
import os
import sys
import traceback

from pytransit.specific_tools import console_tools, logging
from pytransit.globals import gui, cli, root_folder, debugging_enabled
import pytransit.methods


@cli.add_command("help")
def help_command(args=[], kwargs={}):
    from pytransit.generic_tools import misc, command_line
    # 
    # create available subcommand string
    # 
    if True:
        subcommand_prefix      = console_tools.subcommand_prefix
        
        # convert to strings
        subcommands_as_strings = [
            " ".join(each_subcommand)
                for each_subcommand in cli.subcommands.keys()
        ]
        # add color
        subcommands_as_strings_with_color = [
            command_line.color(each_subcommand_string+" ", foreground="bright_green")
                for each_subcommand_string in subcommands_as_strings
        ]
        # add prefix
        subcommands_as_strings_with_color = [
            misc.indent(subcommand_prefix + each_subcommand_string)
                for each_subcommand_string in subcommands_as_strings_with_color
        ]
        # combine
        possible_commands = "\n".join(subcommands_as_strings_with_color)
    
    # 
    # user directly asked for help
    # 
    if len(args) == 0:
        logging.error(
            f"""
                The available subcommands are:\n{possible_commands}
                
                Run a subcommand with no arguments to get subcommand-specific usage/examples
            """.replace("\n                ","\n"),
            no_traceback=True,
        )
    # 
    # user failed at something, and we are trying to help
    # 
    else:
        length_of_longest_subcommand = max(len(each) for each in cli.subcommands.keys())
        given_command                = " ".join(args[:length_of_longest_subcommand])
        given_command_formatted      = subcommand_prefix + command_line.color(" ".join(args[:length_of_longest_subcommand])+" ", foreground="bright_red")
        closest_match, *_            = misc.levenshtein_distance_sort(word=given_command, other_words=subcommands_as_strings)
        closest_match_formatted      = subcommand_prefix + command_line.color(closest_match, foreground="bright_yellow")
        logging.error(
            f"""
                I got this subcommand: {given_command_formatted}
                Maybe you meant:       {closest_match_formatted}
                
                Here are all the available subcommands:\n{possible_commands}
            """.replace("\n                ","\n"),
            no_traceback=True,
        )

# this wrapper exist to help with test cases. Might be good to change the test cases to elimate the need for it
def main(*args, **kwargs):
    # 
    # Check python version
    # 
    if sys.version_info[0] < 3:
        print("TRANSIT v3.0+ requires python3.6+. To use with python2, please install TRANSIT version < 3.0")
        sys.exit(0)
    
    # 
    # parse global args
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
            help_command(args, kwargs)
            sys.exit(0)

    # 
    # GUI Mode
    # 
    if not (args or kwargs):
        # Tried GUI but no wxPython
        if not console_tools.check_if_has_wx():
            print("Please install wxPython to run in GUI Mode. (pip install wxPython)")
            print("")
            print("To run in Console Mode, try 'transit <method>' with one of the following methods:")
            help_command(args, kwargs)
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
        import matplotlib
        matplotlib.use("Agg")
        
        # 
        # try to match longest sequence of arguments with a subcommand
        # 
        length_of_longest_subcommand = max(len(each) for each in cli.subcommands.keys())
        for each_length in reversed(range(1, length_of_longest_subcommand+1)):
            subcommand, subcommand_args = tuple(args[:each_length]), args[each_length:]
            if subcommand in cli.subcommands:
                # if the subcommand exists, run it
                cli.subcommands[subcommand](subcommand_args, kwargs)
                exit(0)
        
        # 
        # runtime only gets here if no subcommands were matched
        # 
        help_command(args, kwargs)

def run_main():
    from pytransit.specific_tools.console_tools import clean_args
    (args, kwargs) = clean_args(sys.argv[1:])
    main(*args, **kwargs)

if __name__ == "__main__":
    run_main()