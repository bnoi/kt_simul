"""
This packages contains the kinetochore dynamics simulation

modules provided:
    - core : main simulation
    - gui: a subpackage with a GUI of the simulation
    - io: input/output functions
"""

__all__ = ["core", "gui", "io"]
__version__ = "1.0"

import logging
import os
import sys

from .utils.color_system import color

def in_ipython():
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True

if in_ipython():
    logformat = '%(asctime)s' + ':'
    logformat += '%(levelname)s' + ':'
    logformat += '%(name)s' + ':'
    # logformat += '%(funcName)s' + ': '
    logformat += ' %(message)s'
else:
    logformat = color('%(asctime)s', 'BLUE') + ':'
    logformat += color('%(levelname)s', 'RED') + ':'
    logformat += color('%(name)s', 'YELLOW') + ':'
    # logformat += color('%(funcName)s', 'GREEN') + ': '
    logformat += color(' %(message)s', 'ENDC')

thisdir = os.path.abspath(os.path.dirname(__file__))
pkgdir = os.path.dirname(thisdir)
samplesdir = os.path.join(pkgdir, 'samples')

logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(logformat, "%Y-%m-%d %H:%M:%S")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
