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
    @staticmethod
    def loess_plot(insertion_counts, ta_sites, output_path=None):
        with gui_tools.nice_error_log:
            import numpy
            import matplotlib
            import matplotlib.pyplot as plt
            from pytransit.specific_tools import stat_tools
            from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
            from pytransit.specific_tools.tnseq_tools import Wig
            from statsmodels.nonparametric.smoothers_lowess import lowess
            
            #
            # get read_counts and positions
            # 
            window = 100
            h      = 20_000
            X, Y = Method.bucketize1(ta_sites=ta_sites, insertion_counts=insertion_counts, window=window)
            smooth      = lambda y_values: stat_tools.smooth_for_loess(X=X, Y=y_values, h=h)
            correct     = lambda y_values: stat_tools.loess_correction2(X=X, Y=y_values, h=h, window=window, should_bucketize=True, should_smooth=False)
            
            Y_smoothed = smooth(Y)
            Y_corrected = correct(Y)
            Y_corrected_and_smoothed = correct(Y_smoothed)
            
            # try plotting with plotly (easier to debug b/c can zoom in/out and hide/show lines)
            if not gui.is_active and misc.go:
                go = misc.go
                fig = go.Figure(
                    layout=misc.go.Layout(
                        title='Multi-Line and Scatterplot Chart',
                        xaxis=dict(title='Genomic Position'),
                        yaxis=dict(title='Mean insertion count over bp windows', type='log'), # log VS linear
                    ),
                )
                fig.add_trace(misc.go.Scatter(x=X, y=Y          , mode='markers', name='Y'),)
                
                group_size = 3
                h_values = [1000,1500,2000,3000,5_000,10_000,20_000,30_000,50_000,70_000,100_000]
                for h in h_values:
                    Y_smoothed = smooth(Y)
                    Y_corrected = correct(Y)
                    Y_corrected_and_smoothed = correct(Y_smoothed)
                    
                    fig.add_trace(misc.go.Scatter(x=X, y=Y_smoothed , mode='lines', name='Y_smoothed'),)
                    fig.add_trace(misc.go.Scatter(x=X, y=Y_corrected, mode='markers', name='Y_corrected', marker=dict(color='#20da9c')),)
                    fig.add_trace(misc.go.Scatter(x=X, y=Y_corrected_and_smoothed, mode='lines', name='Y_corrected_and_smoothed'),)
                
                for index,each in enumerate(fig.data):
                    fig.data[index].visible = False
                fig.data[0].visible = True
                fig.data[1].visible = True
                fig.data[2].visible = True
                fig.data[3].visible = True
                
                # Create and add slider
                steps = []
                for h, index in zip(h_values, range((len(fig.data)-1)//group_size)):
                    actual_index = (index*group_size)+1
                    step = dict(
                        method="update",
                        args=[
                            {"visible": [False] * len(fig.data)},
                            {"title": f"h={h}"},
                        ],  # layout attribute
                    )
                    step["args"][0]["visible"][actual_index] = True  # Toggle i'th trace to "visible"
                    step["args"][0]["visible"][actual_index+1] = True  # Toggle i'th trace to "visible"
                    step["args"][0]["visible"][actual_index+2] = True  # Toggle i'th trace to "visible"
                    steps.append(step)

                sliders = [
                    dict(active=10, currentvalue={"prefix": "h: "}, pad={"t": 50}, steps=steps)
                ]
                
                fig.update_layout(sliders=sliders)
                
                fig.show()
                # misc.go.Figure(
                #     data=[
                #         misc.go.Scatter(x=X, y=Y          , mode='markers', name='Y'),
                #         misc.go.Scatter(x=X, y=Y_smoothed , mode='lines', name='Y_smoothed'),
                #         misc.go.Scatter(x=X, y=Y_corrected, mode='markers', name='Y_corrected', marker=dict(color='#20da9c')),
                #         misc.go.Scatter(x=X, y=Y_corrected_and_smoothed, mode='lines', name='Y_corrected_and_smoothed'),
                #         # misc.go.Scatter(x=X, y=Y_anti_corrected          , mode='markers', name='Y_anti_corrected', marker=dict(color='salmon')),
                #         # misc.go.Scatter(x=X, y=Y_anti_corrected_and_smoothed, mode='lines', name='Y_anti_corrected_and_smoothed'),
                #         # misc.go.Scatter(x=X, y=Y, mode='markers', name='Scatter Points'),
                #     ],
                #     layout=misc.go.Layout(
                #         title='Multi-Line and Scatterplot Chart',
                #         xaxis=dict(title='Genomic Position'),
                #         yaxis=dict(title='Mean insertion count over bp windows', type='log'), # log VS linear
                #     ),
                # ).show()
                
            # fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2)
            # fig, ((ax1, ax3)) = plt.subplots(1, 2)
            fig, ((ax1)) = plt.subplots(1, 1)
            fig.set_figheight(5)
            fig.set_figwidth(10)
            
            ax1.plot(X,Y, "y+")
            ax1.plot(X,Y_smoothed, "b-", label='Before Correction',color="blue")
            ax1.plot(X,Y_corrected, "g-", label='After Correction',color="green")
            ax1.set_xlabel("Genomic Position (bp)")
            ax1.set_ylabel("(log-scale) mean insertion count over %sbp windows" % window)
            ax1.set_yscale("log")
            ax1.set_title("Before Correction")
            plt.legend()
            
            # ax2.plot(X,Y, "g+")
            # ax2.plot(X,fit_no_correction, "b-")
            # ax2.set_xlabel("Genomic Position (bp)")
            # ax2.set_ylabel("mean insertion count over %sbp windows" % window)
            # ax2.set_title("Fit Before Correction")
            
            # ax3.plot(X,Y_corrected, "y+")
            # ax3.plot(X,Y_corrected, "b-")
            # ax3.set_xlabel("Genomic Position (bp)")
            # ax3.set_ylabel("mean insertion count over %sbp windows" % window)
            # ax3.set_yscale("log")
            # ax3.set_title("After LOESS Correction")

            # ax4.plot(X,Y, "y+")
            # ax4.plot(X,fit, "b-")
            # ax4.set_xlabel("Genomic Position (bp)")
            # ax4.set_ylabel("mean insertion count over %sbp windows" % window)
            # ax4.set_title("LOESS Fit")
            
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
    
    # helper
    @staticmethod
    def bucketize1(ta_sites, insertion_counts, window=1000):
        number_of_windows = int(len(ta_sites) / window) + 1  # python3 requires explicit rounding to int
        x_w = numpy.zeros(number_of_windows)
        y_w = numpy.zeros(number_of_windows)
        for window_index in range(number_of_windows):
            x_w[window_index] = window * window_index
            y_w[window_index] = sum(insertion_counts[window * window_index : window * (window_index + 1)])
        return x_w, y_w
        
    @staticmethod
    def bucketize2(x, y, window=1000):
        rounded_x = (x/win).astype(int)*win
        buckets = {}
        for each_bucket, each_count in zip(rounded_x, y):
            if each_bucket not in buckets:
                buckets[each_bucket] = []
            buckets[each_bucket].append(each_count)
        sorted_keys = tuple(sorted(buckets.keys()))
        Y = numpy.array(tuple(numpy.mean(buckets[key]) for key in sorted_keys))
        X = numpy.array(tuple(sorted_keys))
        return X, Y