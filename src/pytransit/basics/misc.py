def open_url(url):
    import sys
    import os
    import subprocess
    
    output = ""
    error = ""
    try:
        if sys.platform.startswith("darwin"): # OSX
            args = ["open", url]
            output, error = subprocess.Popen(
                args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
        elif os.name == "nt": # Windows
            os.startfile(url)
        elif os.name == "posix": # Linux
            args = ["xdg-open", url]
            output, error = subprocess.Popen(
                args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if "not found" in error:
                args = ["exo-open", url]
                output, error = subprocess.Popen(
                    args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()

    except Exception as error:
        raise Exception(f'''Error opening {url} \n{error}''')