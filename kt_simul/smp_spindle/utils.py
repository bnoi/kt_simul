# -*- coding: utf-8 -*-

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function


def indexer(name, index=None, expo=None):
    '''Small utility to append an index  and
    an exponent to the string `name`
    '''
    if index is not None:
        name = name + '_' + index
    if expo is not None:
        name = name + '^' + expo
    return name
