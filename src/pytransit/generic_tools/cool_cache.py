# from: https://raw.githubusercontent.com/andrewgazelka/smart-cache/master/smart_cache/__init__.py
# (MIT License on PyPi)
# has been modified to use super_hash and work on python3.8

import pickle
import queue
import threading
from os import path
from threading import Thread

# from super_hash import super_hash, hash_file
from super_hash import super_hash, hash_file

class Object:
    pass

class NotGiven:
    pass

settings = Object()
settings.default_folder = "cache.ignore/"
settings.worker_que_size = 1000

class CacheData:
    calculated = False
    cache_file_name = ""
    deep_hash = ""
    cache = {}

# since we only care about latest
worker_que = None
def cache(folder=NotGiven, depends_on=lambda:None, watch_attributes=[], watch_filepaths=lambda *args, **kwargs:[], custom_hasher=None, bust=False):
    global worker_que
    
    if folder == NotGiven:
        folder = settings.default_folder
    
    # save in ram
    if folder is None:
        def decorator_name(input_func):
            in_memory_cache = {}
            def wrapper(*args, **kwargs):
                hashed_args = list(args)
                # if watching attributes on self, replace first arg
                if callable(watch_attributes):
                    self = hashed_args[0]
                    hashed_args[0] = watch_attributes(self)
                    
                elif len(watch_attributes) > 0:
                    self = hashed_args[0]
                    attributes = {}
                    for each_attribute in watch_attributes:
                        if hasattr(self, each_attribute):
                            attributes[each_attribute] = getattr(self, each_attribute)
                    hashed_args[0] = attributes
                
                # 
                # filepath hashes
                # 
                filepaths_to_watch = watch_filepaths(*args, **kwargs)
                file_hashes = tuple(hash_file(each) for each in filepaths_to_watch)
                
                # 
                # custom_hasher
                # 
                if callable(custom_hasher):
                    hashed_args = custom_hasher(*args, **kwargs)
                    kwargs = None # need to exclude kwargs when custom_hasher is present
                
                # check if this arg combination has been used already
                arg_hash = super_hash((hashed_args, kwargs, depends_on(), file_hashes))
                if arg_hash in in_memory_cache:
                    return in_memory_cache[arg_hash]
                # if args not in cache, run the function
                else:
                    result = input_func(*args, **kwargs)
                    in_memory_cache[arg_hash] = result # save the output for next time
                    return result
            return wrapper
        return decorator_name
    
    # save in cold storage
    else:
        import file_system_py as FS
        if worker_que is None:
            worker_que = queue.Queue(maxsize=settings.worker_que_size)
            thread = Thread(target=worker)
            thread.start()
        def real_decorator(input_func):
            data = CacheData()  # because we need a reference not a value or compile error
            function_id = super_hash(input_func)
            data.cache_file_name = f'cache.ignore/{function_id}.pickle'
            if bust:
                FS.remove(data.cache_file_name)
            def wrapper(*args, **kwargs):
                hashed_args = list(args)
                # if watching attributes on self, replace first arg
                if callable(watch_attributes):
                    self = hashed_args[0]
                    hashed_args[0] = watch_attributes(self)
                    
                elif len(watch_attributes) > 0:
                    self = hashed_args[0]
                    attributes = {}
                    for each_attribute in watch_attributes:
                        if hasattr(self, each_attribute):
                            attributes[each_attribute] = getattr(self, each_attribute)
                    hashed_args[0] = attributes
                
                # load cached values for this function
                if not data.calculated:
                    data.deep_hash = function_id
                    if path.exists(data.cache_file_name):
                        try:
                            with open(data.cache_file_name, 'rb') as cache_file:
                                func_hash, cache_temp = pickle.load(cache_file)
                                if func_hash == data.deep_hash:
                                    data.cache = cache_temp
                        except Exception as error:
                            # auto remove corrupted files
                            FS.remove(data.cache_file_name)
                    data.calculated = True
                
                # 
                # filepath hashes
                # 
                filepaths_to_watch = watch_filepaths(*args, **kwargs)
                file_hashes = tuple(hash_file(each) for each in filepaths_to_watch)
                
                # 
                # custom_hasher
                # 
                if callable(custom_hasher):
                    hashed_args = custom_hasher(*args, **kwargs)
                    kwargs = None # need to exclude kwargs when custom_hasher is present
                
                # check if this arg combination has been used already
                arg_hash = super_hash((hashed_args, kwargs, depends_on(), file_hashes))
                if arg_hash in data.cache:
                    return data.cache[arg_hash]
                # if args not in cache, run the function
                else:
                    result = input_func(*args, **kwargs)
                    data.cache[arg_hash] = result # save the output for next time
                    data_to_push = CacheData()
                    data_to_push.calculated      = deepcopy(data.calculated)
                    data_to_push.cache_file_name = deepcopy(data.cache_file_name)
                    data_to_push.deep_hash       = deepcopy(data.deep_hash)
                    data_to_push.cache           = deepcopy(data.cache)
                    worker_que.put(data_to_push, block=False) # use a different process for saving to disk to prevent slowdown
                    return result
            return wrapper
        return real_decorator

def worker():
    global worker_que
    while threading.main_thread().is_alive():
        try:
            if worker_que is not None:
                data = worker_que.get(timeout=0.1) # 0.1 second. Allows for checking if the main thread is alive
                while not worker_que.empty(): # so we only write the latest value
                    data = worker_que.get(block=False)
                FS.clear_a_path_for(data.cache_file_name, overwrite=True)
                with open(data.cache_file_name, 'wb') as cache_file:
                    pickle.dump((data.deep_hash, data.cache), cache_file, protocol=4)
                worker_que.task_done()
        except queue.Empty:
            continue