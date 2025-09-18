import os
import sys
import time
import webbrowser
from subprocess import run, Popen
from dirsync import sync

if __name__ == "__main__":
    user_home = os.path.expanduser('~')
    user_nice = os.path.join(user_home, ".nice", "nice_4.0")
    if getattr(sys, 'frozen', False):
        src_nice = os.path.dirname(os.path.dirname(sys.executable))
    else:
        src_nice = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(user_nice):
        os.makedirs(user_nice)
        print("created", user_nice)
    print("syncing", src_nice, "->", user_nice)
    sync(src_nice, user_nice, 'sync', exclude=('^.git', '^build','^db.sqlite3', '.*__pycache__.*'))
    venv_dir = os.path.join(user_nice, "venv")
    if not os.path.exists(venv_dir):
        print("creating virtual python environment", venv_dir)
        run(["python3", "-m", "venv", "venv"], cwd=user_nice)
    print("updating virtual python environment", venv_dir)
    run(["bin/python3", "-m", "pip", "install", "--upgrade", "pip"], cwd=venv_dir)
    print("installing dependencies into virtual python environment", venv_dir)
    run(["bin/pip", "install", "Django==5.2.5"], cwd=venv_dir)
    print("updating database")
    run(["venv/bin/python3", "manage.py", "makemigrations"], cwd=user_nice)
    run(["venv/bin/python3", "manage.py", "migrate"],        cwd=user_nice)
    print("starting server")
    server = Popen(["venv/bin/python3", "manage.py", "runserver", "--noreload"], cwd=user_nice)
    print("opening URL", "http://localhost:8000")
    webbrowser.open_new("http://localhost:8000")
    server.communicate() # block until server is killed (ctrl-C)
