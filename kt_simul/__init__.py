"""
This packages contains the kinetochore dynamics simulation

modules provided:
    - core : main simulation
    - gui: a subpackage with a GUI of the simulation
    - io: input/output functions
"""

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

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

log = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(logformat, "%Y-%m-%d %H:%M:%S")
handler.setFormatter(formatter)
log.addHandler(handler)
log.setLevel(logging.DEBUG)
log.propagate = False

log = logging.getLogger(__name__)

try:
    import tqdm
except:
    log.warning("For progress bar you need to pip install tqdm")
