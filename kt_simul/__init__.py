"""
This packages contains the kinetochore dynamics simulation

modules provided:
    - core : main simulation
    - gui: a subpackage with a GUI of the simulation
"""
__all__ = ["core", "gui"]

import logging
import os
import sys
import numpy

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

import pyximport

# Compile Cython file
if os.name == 'nt':
    if os.environ.has_key('CPATH'):
        os.environ['CPATH'] = os.environ['CPATH'] + numpy.get_include()
    else:
        os.environ['CPATH'] = numpy.get_include()

    # XXX: we're assuming that MinGW is installed in C:\MinGW (default)
    if os.environ.has_key('PATH'):
        os.environ['PATH'] = os.environ['PATH'] + ';C:\MinGW\bin'
    else:
        os.environ['PATH'] = 'C:\MinGW\bin'

    mingw_setup_args = { 'options': { 'build_ext': { 'compiler': 'mingw32' } } }
    pyximport.install(setup_args=mingw_setup_args)

elif os.name == 'posix':
    if os.environ.has_key('CFLAGS'):
        os.environ['CFLAGS'] = os.environ['CFLAGS'] + ' -I' + numpy.get_include()
    else:
        os.environ['CFLAGS'] = ' -I' + numpy.get_include()

    pyximport.install()
