#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This packages contains the kinetochore dynamics simulation

modules provided:
    - core : main simulation
    - gui: a subpackage with a GUI of the simulation
    - analysis : Data processing, analysis, batch scripts
"""
__all__ = ["core", "analysis", "gui"]

import logging
import os
import sys
import numpy

from kt_simul.utils.color import color

# Setup logging

logformat = color('%(asctime)s', 'BLUE') + ':'
logformat += color('%(levelname)s', 'RED') + ':'
logformat += color('%(name)s', 'YELLOW') + ':'
logformat += color('%(message)s', 'ENDC')
logging.basicConfig(level=logging.DEBUG,
                    format=logformat,
                    datefmt='%Y-%m-%d %H:%M:%S')

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
