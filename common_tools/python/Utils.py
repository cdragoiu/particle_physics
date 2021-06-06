import os

# implementation of the 'which' command ------------------------------------------------------------
def Which(file):
    for path in os.environ['PATH'].split(':'):
        try:
            if file in os.listdir(path):
                return path + '/' + file
        except OSError:
            continue
