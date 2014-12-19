# -*- coding: utf-8 -*-
"""
Handler for the parameters files
"""

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import logging
import io

import numpy as np
import pandas as pd

# Those strings should be respected in the xml file
SPRING_UNIT = u'pN/µm'
DRAG_UNIT = u'pN.s/µm'
LENGTH_UNIT = u'µm'
FREQ_UNIT = u'Hz'
FORCE_UNIT = u'pN'
SPEED_UNIT = u'µm/s'

__all__ = ["ParamTree"]

log = logging.getLogger(__name__)


class ParamTree(object):
    """
    This class defines the container for the simulation parameters.
    It wraps an ElementTree instance whose elements contains the
    name, value (as a string), description and unit of each parameter and a
    dictionnary that is used during the simulation.

    The `value` attribute of the tree is not modified by adimentionalization,
    whereas the value in the dictionnary is changed.
    """

    def __init__(self, df=None, jsonfile=None, adimentionalized=True):
        """
        """
        if jsonfile:
            self.params = pd.read_json(jsonfile, orient='records')
        elif isinstance(df, pd.DataFrame):
            self.params = df

        self.adimentionalized = adimentionalized
        self.build_tree()

    def build_tree(self):
        """
        """
        self.absolute_dic = self.params.loc[:, ['name', 'value']]
        self.absolute_dic = self.absolute_dic.set_index('name').to_dict()['value']
        self.relative_dic = self.absolute_dic.copy()

        if self.adimentionalized:
            self.adimentionalize()

    def adimentionalize(self):
        """
        This function scales everything taking dt as unit time, Vk as
        unit speed, and Fk as unit force. It relies on a correct
        definition of the units of the elements of the param tree,
        thus a correct spelling in the json file, so please beware.
        """
        Vk = self.absolute_dic["Vk"]
        Fk = self.absolute_dic["Fk"]
        dt = self.absolute_dic["dt"]

        for i, p in self.params.iterrows():
            key = p['name']
            value = p['value']
            if p['unit'] == SPRING_UNIT:
                value /= Fk
            elif p['unit'] == DRAG_UNIT:
                value *= Vk / Fk
            elif p['unit'] == SPEED_UNIT:
                value /= Vk
            elif p['unit'] == FREQ_UNIT:
                value *= dt
            elif p['unit'] == FORCE_UNIT:
                value /= Fk
            self.relative_dic[key] = value

    def save(self, path):
        """Save as json file
        """
        f = open(path, "w")
        f.write(self.params.to_json(orient='records'))
        f.close()

    def __setitem__(self, key, item):
        """
        """
        params = self.params.set_index('name')
        params.loc[key, 'value'] = item
        self.params = params.reset_index()

        self.build_tree()
