import sys
import os
import time
import ntpath
import math
import random
import datetime
import collections
import heapq

import numpy

from pytransit.specific_tools import  gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import csv, misc
from pytransit.specific_tools.transit_tools import wx, basename
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

@misc.singleton
class Method:
    name = "LOESS Plot"
    usage_string = f"""
        Usage:
            {console_tools.subcommand_prefix} loess <wig_file> <where_to_save.png>
    """.replace("\n        ", "\n")
    
    # 
    # main method
    # 
    # currently generates plots for only first sample (assumes K=1)
    #
    @staticmethod
    def loess_plot(insertion_counts, ta_sites, output_path=None):
        with gui_tools.nice_error_log:
            import numpy
            import matplotlib
            import matplotlib.pyplot as plt
            from pytransit.specific_tools import stat_tools
            from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
            from pytransit.specific_tools.tnseq_tools import Wig
            import pytransit.specific_tools.stat_tools
                
            #
            # get read_counts and positions
            #
            data, position = insertion_counts,ta_sites
            if isinstance(data, (list, tuple)):
                data = numpy.array(data)
            
            if data.ndim==1: data = data[:,None] # vec->arr, in case there is only 1 sample
            data = data.transpose()
            (K,N) = data.shape
            window = 100

            # assume K=1 (show plots for first sample only)
            size = int(len(position)/window) + 1 # python3 requires explicit rounding to int
            x_w = numpy.zeros(size)
            y_w = numpy.zeros(size)
            for i in range(size):
                x_w[i] = window*i
                y_w[i] = sum(data[0][window*i:window*(i+1)])

            fig,((ax1,ax2)) = plt.subplots(1,2)
            fig.set_figheight(5)
            fig.set_figwidth(11)
            
            y_smooth = pytransit.specific_tools.stat_tools.smooth_for_loess(x_w, y_w, h=10000)
            ax1.plot(x_w, y_w, "g+")
            ax1.plot(x_w, y_smooth, "b-")
            ax1.set_xlabel("Genomic Position (TA sites)")
            ax1.set_ylabel("Reads per 100 insertion sites")
            
            # assume K=1
            y2 = pytransit.specific_tools.stat_tools.loess_correction(position,data[0],h=10000)
            size = int(len(position)/window) + 1 # python3 requires explicit rounding to int
            y2_w = numpy.zeros(size)
            for i in range(size):
                y2_w[i] = sum(y2[window*i:window*(i+1)])

            y2_smooth = pytransit.specific_tools.stat_tools.smooth_for_loess(x_w, y2_w, h=10000)
            ax2.plot(x_w, y2_w, "g+")
            ax2.plot(x_w, y2_smooth, "b-")
            ax2.set_xlabel("Genomic Position (TA sites)")
            ax2.set_ylabel("Reads per 100 insertion sites")
            ax2.set_title("After LOESS correction")
            

            if output_path == None:
                plt.show()
            else:
                if not output_path.endswith(".png"):
                    output_path = output_path+".png"
                print(f"saving to: {output_path}")
                plt.savefig(output_path)

    # 
    # CLI
    #
    @staticmethod 
    @cli.add_command("loess")
    def loess_cli(args, kwargs):
        from pytransit.specific_tools.tnseq_tools import Wig
        
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=2)
        
        wig_file_path   = args[0]
        output_path     = args[1]
        wig = Wig.read(wig_file_path)
        
        Method.loess_plot(
            insertion_counts=wig.insertion_counts,
            ta_sites=wig.ta_sites,
            output_path=output_path,
        )
            
    # 
    # GUI
    # 
    @gui.add_wig_area_dropdown_option(name=name)
    def click_show_loess(event):
        import numpy
        from pytransit.specific_tools.tnseq_tools import Wig
        
        # 
        # get selection
        # 
        wig_objects = gui.selected_samples  or  gui.samples
        
        #
        # get read_counts and positions
        # 
        insertion_counts_matrix, ta_sites = Wig.selected_as_gathered_data(wig_objects)
        Method.loess_plot(
            insertion_counts=numpy.mean(insertion_counts_matrix,axis=0),
            ta_sites=ta_sites,
        )