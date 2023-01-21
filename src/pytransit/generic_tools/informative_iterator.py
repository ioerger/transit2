# heavily adpted from: https://raw.githubusercontent.com/PaoloLRinaldi/progress_bar_python/master/perc.py
from datetime import datetime, timedelta
from statistics import mean, stdev
import time
import sys
import math

from pytransit.generic_tools.lazy_dict import LazyDict as Map

# nested indentation support
try:
    from blissful_basics import print as bliss_print
except Exception as error:
    # create a placeholder/stand-in
    def bliss_print(*args, **kwargs):
        print(*args, **kwargs)
    class Indent:
        size = 0
        string = " "
    bliss_print.indent = Indent()

# GUI progress support in notebooks
try:
    from IPython.display import display, HTML, clear_output
    from io import StringIO
    ipython_exists = False
    
    try:
        from IPython.display import get_ipython
        ipython_exists = 'IPKernelApp' in get_ipython().config
        ipython_exists = True
    except Exception as error:
        pass
    
    try:
        from google.colab import output
        ipython_exists = True
    except Exception as error:
        pass
    
except ImportError:
    ipython_exists = False
except AttributeError:
    ipython_exists = False
    
class NotGiven: pass

def subsequence_replace(a_list, sequence, replacement):
    que = []
    new_list = []
    for each in a_list:
        que.append(each)
        while len(que) > len(sequence):
            new_list.append(que.pop(0))
        if que == sequence:
            que.clear()
            for each in replacement:
                new_list.append(each)
    new_list += que # add any remaining elements
    return new_list    

class ProgressBar:
    """
    from informative_iterator import ProgressBar, time
    for progress, each in ProgressBar(10000):
        time.sleep(0.01)
    """
    
    layout = [ 'title', 'bar', 'percent', 'spacer', 'fraction', 'spacer', 'remaining_time', 'spacer', 'end_time', 'spacer', 'duration', 'spacer', ]
    minimal_layout = [ 'title', 'bar', 'spacer', 'end_time', 'spacer', ]
    spacer = " | "
    minmal = False
    inline = True
    disable_logging = False
    progress_bar_size = 35
    seconds_per_print = 0.1
    percent_per_print = 2
    lookback_size = 100
    time_format = "%H:%M:%S"
    long_time_format = "%D %H:%M:%S"
    
    @classmethod
    def configure(this_class, **config):
        for each_key, each_value in config.items():
            setattr(this_class, each_key, each_value)
    
    def __init__(self, iterator, *, title=None, iterations=None, layout=None, disable_logging=NotGiven, minimal=NotGiven, inline=NotGiven, progress_bar_size=NotGiven, seconds_per_print=NotGiven, percent_per_print=NotGiven, lookback_size=NotGiven):
        original_generator = range(int(iterator)) if isinstance(iterator, (int, float)) else iterator
        self.title = title
        
        # inherit unspecified options from class object
        for each_option in [ "disable_logging", "minimal", "inline", "progress_bar_size", "seconds_per_print", "percent_per_print", "lookback_size" ]:
            arg_value = eval(each_option, locals())
            # default to the class value if not given
            if arg_value == NotGiven:
                actual_value = getattr(ProgressBar, each_option, None)
            # otherwise use the given value
            else:
                actual_value = arg_value
            # set the object's value
            setattr(self, each_option, actual_value)
        
        # initilize misc values
        self.past_indicies     = []
        self.start_time        = datetime.now()
        self.next_percent_mark = self.percent_per_print
        self.prev_time         = -math.inf
        self.times             = [time.time()]
        self.opacity = 0.7
        self.colors = Map(
            progress="#9b68ab",
        )
        self.progress_data     = Map(
            index=0,
            percent=0,
            updated=True,
            time=self.times[0],
            total_iterations=(len(original_generator) if iterations is None else iterations),
            deviation=None,
            expected_number_of_updates_needed=None,
            pretext="",
            text="",
        )
        # setup print
        if self.disable_logging:
            self.print = lambda *args, **kwargs: None
        elif not ipython_exists:
            self.print = print
        else:
            # remove the progress bar and percent
            layout = list(self.layout)
            layout = subsequence_replace(layout, [ 'spacer', 'bar'    , 'spacer' ], ['spacer'])
            layout = subsequence_replace(layout, [ 'spacer', 'bar'    ,          ], ['spacer'])
            layout = subsequence_replace(layout, [           'bar'    , 'spacer' ], ['spacer'])
            layout = subsequence_replace(layout, [           'bar'               ], [])
            layout = subsequence_replace(layout, [ 'spacer', 'percent', 'spacer' ], ['spacer'])
            layout = subsequence_replace(layout, [ 'spacer', 'percent',          ], ['spacer'])
            layout = subsequence_replace(layout, [           'percent', 'spacer' ], ['spacer'])
            layout = subsequence_replace(layout, [           'percent',          ], [])
            layout = subsequence_replace(layout, [ 'spacer', 'spacer'            ], ['spacer'])
            self.layout = layout
            self.string_buffer = ""
            self.html_created = False
            self.should_flush = False
            def ipython_print(*args, **kwargs):
                # get the string value
                string_stream = StringIO()
                print(*args, **kwargs, file=string_stream)
                output_str = string_stream.getvalue()
                string_stream.close()
                self.string_buffer += output_str
                # escape html just encase
                self.string_buffer = self.string_buffer.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('‘', "'").replace('"', "&quot;").replace("\n", "<br>")
                
                # clear output whenever newline is created
                if self.should_flush:
                    self.should_flush = False
                    if not self.html_created:
                        self.html_created = True
                        clear_output(wait=True)
                        display(HTML(f'''
                            <div id="progressContainer" style="left: 0; top: 0; width: 95%; background: transparent; color: white; position: sticky; padding: 1.2rem 0.4rem; box-sizing: border-box;">
                                <div style="position: relative; background: #46505a; height: 2.7rem; width: 100%; border-radius: 10rem; border: transparent solid 0.34rem; box-sizing: border-box; opacity: {self.opacity};">
                                    <!-- color bar -->
                                    <div id="progressBar" style="height: 100%; background: {self.colors.progress}; border-radius: 10rem; transition: all 0.5s ease-in-out 0s;"></div>
                                    <!-- percentage text -->
                                    <div style="height: 100%; width: 100%; position: absolute; top: 0; left: 0; display: flex; flex-direction: column; align-items: center; align-content: center; justify-items: center;  justify-content: center;">
                                        <span id="progressPercent">
                                            {self.progress_data.percent:0.2f}%
                                        </span>
                                    </div>
                                </div>
                                <div style="width: 100%; height: 1rem;">
                                </div>
                                <div style="position: relative; height: fit-content; width: 100%; box-sizing: border-box; display: flex; flex-direction: column; align-items: center; align-content: center; justify-items: center;  justify-content: center;">
                                    <div style="position: relative;background: #46505a;height: fit-content;width: fit-content; min-width: 50%; border-radius: 1.2rem;border: transparent solid 0.5rem;box-sizing: border-box;display: flex;flex-direction: column;align-items: center;align-content: center;justify-items: center;justify-content: center;padding: 1rem; padding-top: 0;">
                                        <code id="progressText" style="white-space: pre; color: whitesmoke;" >{self.string_buffer}</code>
                                        <code id="customProgressText" style="white-space: pre; color: whitesmoke; padding-top: 1rem;" >{self.progress_data.text}</code>
                                    </div>
                                </div>
                                <style>
                                </style>
                                <script>
                                    var progressContainer = document.getElementById("progressContainer")
                                    var progressFooter = document.getElementById("progressFooter")
                                    var outputArea = document.getElementById("output-area")
                                    var outputHeader = document.getElementById("output-header")
                                    var outputBody = document.getElementById("output-body")
                                    outputHeader.appendChild(progressContainer)
                                    
                                    var scrollBoxSize = 27
                                    // outputArea
                                    outputArea.style.maxHeight = `${"{scrollBoxSize+10}"}rem`
                                    outputArea.style.minHeight = `${"{scrollBoxSize+10}"}rem`

                                    // outputHeader
                                    outputHeader.style.position = "absolute"
                                    outputHeader.style.top = "0"
                                    outputHeader.style.left = "0"
                                    outputHeader.style.width = "100%"

                                    // outputBody
                                    outputBody.style.paddingTop = "10.5rem"
                                    outputBody.style.position = "absolute"
                                    outputBody.style.height = `${"{scrollBoxSize}"}rem`
                                    outputBody.style.overflow = "auto"
                                    outputBody.style.minHeight = "50vh"
                                    outputBody.style.top = "0"
                                    outputBody.style.width = "97%"

                                    outputBody.scrollTo(0, outputBody.scrollHeight)
                                </script>
                            </div>
                        '''))
                    else:
                        from random import random
                        random_id_1 = f"id_{random()}".replace('.','')
                        random_id_2 = f"id_{random()}".replace('.','')
                        display(HTML(f'''
                            <div>
                                <code id="{random_id_1}" style="display: none;" >{self.string_buffer}</code>
                                <code id="{random_id_2}" style="display: none;" >{self.progress_data.text}</code>
                                <script>
                                    var bar = document.getElementById("progressBar")
                                    var percent = document.getElementById("progressPercent")
                                    var text = document.getElementById("progressText")
                                    var customText = document.getElementById("customProgressText")
                                    var textContainer = document.getElementById("{random_id_1}")
                                    var customTextContainer = document.getElementById("{random_id_2}")
                                    bar.style.width = `{self.progress_data.percent}%`
                                    percent.innerHTML = `{self.progress_data.percent:0.2f}%`
                                    // swap out contents (performed this way so that self.string_buffer uses html-escapes instead of javascript-string escapes)
                                    text.innerHTML = textContainer.innerHTML
                                    if (customTextContainer.innerHTML.trim().length > 0) customText.innerHTML = customTextContainer.innerHTML
                                </script>
                            </div>
                        '''))
                        
                    # clear the buffer
                    self.string_buffer = ""
            self.print = ipython_print
            
        # setup layout
        if layout == None and self.minimal:
            self.layout = ProgressBar.minimal_layout
        elif layout == None:
            self.layout = ProgressBar.layout
        else:
            self.layout = layout
        
        def generator_func():
            for self.progress_data.index, each_original in enumerate(original_generator):
                # update
                self.progress_data.time    = time.time()
                self.progress_data.updated = (self.progress_data.time - self.prev_time) > self.seconds_per_print
                self.progress_data.percent = (self.progress_data.index * 10000 // self.progress_data.total_iterations) / 100 # two decimals of accuracy
                
                if self.progress_data.updated:
                    self.prev_time = self.progress_data.time
                    self.times.append(self.progress_data.time)
                    self.past_indicies.append(self.progress_data.index)
                    
                # also printout at each percent marker
                if self.progress_data.percent >= self.next_percent_mark:
                    self.next_percent_mark += self.percent_per_print
                    self.progress_data.updated = True
                
                if self.progress_data.updated:
                    self.total_eslaped_time = 0
                    self.eslaped_time       = 0
                    self.secs_remaining     = math.inf
                    # compute changes (need at least two times to do that)
                    if len(self.times) > 3:
                        self.total_eslaped_time = self.times[-1] - self.times[ 0]
                        self.eslaped_time       = self.times[-1] - self.times[-2]
                        
                        # 
                        # compute ETA as a slight overestimate that is less-of-an-overesitmate over time
                        # 
                        lookback_size = self.lookback_size
                        remaining_number_of_iterations = self.progress_data.total_iterations - self.progress_data.index
                        recent_indicies                = self.past_indicies[-lookback_size:]
                        recent_update_times            = self.times[-lookback_size:]
                        number_of_indicies_processed   =     recent_indicies[-1] -     recent_indicies[0]
                        time_span_of_recent_updates    = recent_update_times[-1] - recent_update_times[0]
                        time_per_update                = time_span_of_recent_updates / (len(recent_update_times)-1)
                        
                        list_of_iterations_per_update = tuple(each-prev for prev, each in zip(recent_indicies[0:-1], recent_indicies[1:]))
                        average_number_of_iterations_per_update = mean(list_of_iterations_per_update)
                        stdev_of_iters_per_update               = stdev(list_of_iterations_per_update) # TODO: the proper way to do this would be with a one sided bell curve
                        partial_deviation                       = (1 - (self.progress_data.percent/100)) * stdev_of_iters_per_update
                        iterations_per_update_lowerbound        = max((average_number_of_iterations_per_update-partial_deviation, min(list_of_iterations_per_update)))
                        soft_lowerbound                         = mean((iterations_per_update_lowerbound, average_number_of_iterations_per_update))
                        expected_number_of_updates_needed       = remaining_number_of_iterations / soft_lowerbound
                        
                        self.secs_remaining = time_per_update * expected_number_of_updates_needed
                    
                    indent = bliss_print.indent.string*bliss_print.indent.size
                    if self.progress_data.pretext:
                        self.print('', end='\r')
                        self.print('                                                                                                                        ', end='\r')
                        bliss_print(self.progress_data.pretext, end='\n')
                    
                    if self.inline:
                        self.print('', end='\r')
                    
                    self.print(indent, end="")
                    
                    # display each thing according to the layout
                    for each in self.layout:
                        getattr(self, f"show_{each}", lambda : None)()
                    
                    if not ipython_exists:
                        if self.progress_data.text:
                            line_1, *other_lines = self.progress_data.text.split('\n')
                            self.print(line_1, end='')
                            if other_lines:
                                nextline_text = f"\n{indent}" + f"\n{indent}".join(other_lines)
                                self.print(nextline_text, end='')
                        
                    if not self.inline:
                        self.print()
                    
                    # convoluted so that it handles both GUI and CLI
                    self.should_flush = True
                    self.print(end="")
                    sys.stdout.flush()
                    
                yield self.progress_data, each_original
                # manual stop if given an infinite generator and a "total_iterations" argument
                if self.progress_data.index+1 >= self.progress_data.total_iterations:
                    break
            
            self.show_done()
            
        self.iterator = iter(generator_func())

    def to_time_string(self, secs):
        secs  = int(round(secs))
        mins  = secs  // 60; secs  = secs  % 60
        hours = mins  // 60; mins  = mins  % 60
        days  = hours // 24; hours = hours % 24
        if days:
            hours = f"{hours}".rjust(2,"0")
            mins  = f"{mins}".rjust(2,"0")
            return f"{days}days, {hours}:{mins}min"
        elif hours:
            mins  = f"{mins}".rjust(2,"0")
            return f"{hours}:{mins}min"
        elif mins:
            secs  = f"{secs}".rjust(2,"0")
            return f"{mins}:{secs}sec"
        else:
            return f"{secs}sec"
            
    def show_spacer(self):
        self.print(self.spacer, end='')
    
    def show_title(self):
        if self.title is not None:
            self.print(self.title, end=' ')
        
    def show_bar(self):
        prog = int((self.progress_data.index / self.progress_data.total_iterations) * self.progress_bar_size)
        self.print('[' + '=' * prog, end='')
        if prog != self.progress_bar_size:
            self.print('>' + '.' * (self.progress_bar_size - prog - 1), end='')
        self.print('] ', end='')
    
    def show_remaining_time(self):
        if self.secs_remaining == math.inf:
            self.print(f'remaining: _______', end='')
        elif self.progress_data.percent != 100:
            self.print(f'remaining: {self.to_time_string(self.secs_remaining)}', end='')
        
    def show_percent(self):
        self.print(f'{self.progress_data.percent:.2f}%'.rjust(6), end='')
    
    def show_duration(self):
        self.print("elapsed: "+self.to_time_string(self.total_eslaped_time), end='')
    
    def show_fraction(self):
        total_str = f"{self.progress_data.total_iterations}"
        self.print(f'{self.progress_data.index}'.rjust(len(total_str))+f'/{self.progress_data.total_iterations}', end='')
    
    # TODO: fix then add this back
    # def show_iteration_time(self):
    #     iterations_per_sec = (self.past_indicies[-1] - self.past_indicies[-2]) / self.eslaped_time
    #     self.print(f'{self.to_time_string(self.eslaped_time)}sec per iter', end='')
    
    def show_start_time(self):
        if self.progress_data.percent != 100:
            self.print(f'started: {self.start_time.strftime("%H:%M:%S")}',  end='')
    
    def show_end_time(self):
        if self.progress_data.percent != 100:
            time_format = self.time_format
            if self.secs_remaining > (86400/2): # more than half a day
                time_format = self.long_time_format
            
            try:
                endtime = self.start_time + timedelta(seconds=self.total_eslaped_time + self.secs_remaining)
                self.print(f'eta: {endtime.strftime(time_format)}',  end='')
            except:
                self.print(f'eta: {"_"*(len(time_format)-3)}',  end='')
    
    def show_done(self):
        if self.inline:
            print("")
        duration = self.to_time_string(time.time() - self.times[0])
        end_time = datetime.now().strftime("%H:%M:%S")
        self.progress_data.percent = 100.0
        self.string_buffer = "" # for ipython
        indent = bliss_print.indent.string * bliss_print.indent.size
        self.print(f'{indent}Done in {duration} at {end_time}')

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.iterator)
    
    def __len__(self):
        return self.progress_data.total_iterations
