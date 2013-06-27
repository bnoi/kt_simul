from __future__ import with_statement
from fabric.api import *
from fabric.colors import *

import os
import socket

# Local paths
project = os.path.abspath(os.path.dirname(__file__))
project = os.path.dirname(os.path.dirname(project))
if socket.gethostname() == "boromir":
    results = "/home/hadim/local/data/"
else:
    results = "/media/thor/data/ktsimu/"

# Loki paths
host = 'hadim@130.120.107.234'
rproject = "/home/hadim/dev/kt_simul/"
rresults = "/home/hadim/local/data/ktsimu/"
rpython = "/home/hadim/local/virtualenvs/ktsimu/bin/python "
rpythonpath = "/home/hadim/dev/kt_simul"

env.hosts = [host, ]

NSIMU = 10000

@task
def push():
    """
    Push code to host
    """
    print red("Push kt_simul to loki")
    with lcd(project):
        local("rsync --progress -a --delete ../kt_simul/ %s:%s" % (host, rproject))

@task
def pull():
    """
    """
    cmd = "rsync --exclude raw --progress -a %s:%s %s" % (host, rresults, results)
    local(cmd)


# Function below are not really used anymore
# Instead I launch python scripts see in scripts/loki/ folder

@task
def launch(nsimu=None, name=""):
    """
    Launch simulations
    """
    if not nsimu:
        nsimu = NSIMU
    push()
    with cd(os.path.join(rproject, "scripts")):
        run("workon ktsimu")
        cmd = rpython + "cluster.py --path %s --nsimu %s --name %s" \
            % (rresults, str(nsimu), name)
        # Allow to run in background
        run("screen -dmS ktsimu " + cmd)


@task
def pool(path):
    """
    Launch PoolEvaluations on simulations
    """
    push()
    with cd(os.path.join(rproject, "scripts")):
        run("workon ktsimu")
        cmd = rpython + "pool.py --path %s" % (path)
        # Allow to run in background
        run(cmd)
        # run("screen -dmS ktsimu " + cmd)


@task
def kill():
    """
    Ugly way : kill all python process
    """
    run("killall -9 python")


@task
def status(keep=False):
    """
    Display the status of all running simulations
    """
    for simu in _listdir(rresults):
        if simu:
            simu_path = os.path.join(rresults, simu)
            files = _listdir(simu_path)
            # If no simu.log then simu is running
            if not "simu.log" in files:
                _show_status(simu_path, keep)


def _listdir(path):
    files = run("ls %s" % path).replace("\n", " ")
    files = files.replace("\t", " ").replace("\r", " ").split(" ")
    return files

def _show_status(path, keep):
    print red("Simulation %s" % path)
    if keep:
        last = " -f "
    else:
        last = " -1 "
    run("tail %s %s" % (last, os.path.join(path, "run.log")))
