from time import time as now
prev = now()
def time_since_prev(*args):
    global prev
    current = now()
    print(*args, "duration",f"{current-prev:0.5f}")
    prev = current