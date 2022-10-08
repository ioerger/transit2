import sys

foreground_colors = {
    "black"          : 30,
    "red"            : 31,
    "green"          : 32,
    "yellow"         : 33,
    "blue"           : 34,
    "magenta"        : 35,
    "cyan"           : 36,
    "white"          : 37,
    "bright_black"   : 90,
    "bright_red"     : 91,
    "bright_green"   : 92,
    "bright_yellow"  : 93,
    "bright_blue"    : 94,
    "bright_magenta" : 95,
    "bright_cyan"    : 96,
    "bright_white"   : 97,
}
background_colors = {
    "black"          : 40,
    "red"            : 41,
    "green"          : 42,
    "yellow"         : 43,
    "blue"           : 44,
    "magenta"        : 45,
    "cyan"           : 46,
    "white"          : 47,
    "bright_black"   : 100,
    "bright_red"     : 101,
    "bright_green"   : 102,
    "bright_yellow"  : 103,
    "bright_blue"    : 104,
    "bright_magenta" : 105,
    "bright_cyan"    : 106,
    "bright_white"   : 107,
}

def color(string, foreground="white", background="black"):
    # if outputing to a file, dont color anything
    if not sys.stdout.isatty():
        return string
        
    if foreground not in foreground_colors:
        raise Exception(f"couldn't find foreground color {foreground}")
    if background not in background_colors:
        raise Exception(f"couldn't find background_color color {background}")
    foreground_number = foreground_colors[foreground]
    background_number = background_colors[background]
    
    foreground_code = f"\u001b[{foreground_number}m"
    background_code = f"\u001b[{background_number}m"
    reset_code = f"\u001b[0m"
    return foreground_code+background_code+str(string)+reset_code
