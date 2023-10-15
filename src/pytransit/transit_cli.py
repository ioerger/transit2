#!/usr/bin/env python3
import os
import sys
import traceback

from pytransit.specific_tools import console_tools
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

def parse_global_flags(args, kwargs):
    # 
    # extract debug value
    # 
    import pytransit.globals as transit_globals
    transit_globals.debugging_enabled = "--debug" in sys.argv
    if transit_globals.debugging_enabled:
        args.remove("--debug")
        sys.argv.remove("--debug")
        kwargs.pop("-debug")
    
    # 
    # special options
    # 
    if len(args) == 0 and len(kwargs) > 1:
        if "help" in kwargs or "h" in kwargs:
            intentional_help_command([], kwargs)
        elif "version" in kwargs or "v" in kwargs:
            version_command([],kwargs)
    
def find_and_run_subcommand(args, kwargs):
    length_of_longest_subcommand = max(len(each) for each in cli.subcommands.keys())
    for each_length in reversed(range(1, length_of_longest_subcommand+1)):
        subcommand, subcommand_args = tuple(args[:each_length]), args[each_length:]
        if subcommand in cli.subcommands:
            # if the subcommand exists, run it
            cli.subcommands[subcommand](subcommand_args, kwargs)
            exit(0)

@cli.add_command("help")
def intentional_help_command(args, kwargs):
    help_command(args, kwargs)
    exit(0)

def help_command(args=[], kwargs={}):
    from pytransit.generic_tools import misc, command_line
    # 
    # create available subcommand string
    # 
    if True:
        # subcommand_prefix      = console_tools.subcommand_prefix
        subcommand_prefix = ""
        
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
            exit_code=0,
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

@cli.add_command("version")
def version_command(args, kwargs):
    import pytransit
    print("Version: {0}".format(pytransit.__version__))
    exit(0)