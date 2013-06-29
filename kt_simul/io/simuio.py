"""
Main module for save and read simulation.
"""

import time
import logging
import StringIO
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from kt_simul.core.spindle_dynamics import KinetoDynamics
from kt_simul.core.simul_spindle import Metaphase
from kt_simul.io.xml_handler import ParamTree, indent, ResultTree

logger = logging.getLogger(__name__)


class SimuIO():
    """
    Parameters
    ----------
    meta_instance : :class:`~kt_simul.core.simul_spindle.Metaphase`
        Should have been already perform a simulation

    """

    def __init__(self, meta_instance=None):

        if meta_instance:
            self.meta = meta_instance
            self.KD = self.meta.KD
            self.paramtree = self.meta.paramtree
            self.measuretree = self.meta.measuretree

    def save(self, simufname, metadata=None, verbose=False):
        """
        Save :class:`~kt_simul.core.simul_spindle.Metaphase` instance to
        HDF5 file. Each of the following element will be stored via
        :class:`pandas.DataFrame` or :class:`pandas.Series` object.

        Parameters
        ----------
        simufname : str
            Filename to save HDF5 file.
        metadata : dict
            Will be stored as a :class:`pandas.Series`. (Not implemented)
        verbose : bool
            Enable verbose mode.

        """

        KD = self.meta.KD
        paramtree = self.meta.paramtree
        measuretree = self.meta.measuretree
        timelapse = self.meta.timelapse
        chromosomes = self.KD.chromosomes

        time_index = pd.MultiIndex.from_arrays([timelapse], names=['t'])

        df_to_save = ['params', 'measures', 'spbs', 'kts', 'plug_sites']

        """
        spbs DataFrame look like that:

                               x
        t   label side
        0.0 spb   A     0.150000
            spbL  B    -0.150000
        0.5 spb   A     0.154235
            spbL  B    -0.154235
        1.0 spb   A     0.158469
            spbL  B    -0.158469
        1.5 spb   A     0.162702
            spbL  B    -0.162702
        2.0 spb   A     0.166934
            spbL  B    -0.166934
        """
        spbR = pd.DataFrame(KD.spbR.traj, columns=['x'], index=time_index)
        spbL = pd.DataFrame(KD.spbL.traj, columns=['x'],index=time_index)

        spbs = pd.concat({('spb', 'A'): spbR, ('spb', 'B'): spbL}, names=['label', 'side'])
        spbs = spbs.reorder_levels([2, 0, 1]).sort()

        """
        kts DataFrame look like that:

                                  x
        t   id label side
        0.0 0  kt    A    -0.077064
                     B     0.022936
            1  kt    A    -0.037428
                     B     0.062572
            2  kt    A    -0.136443
                     B    -0.036443
        0.5 0  kt    A    -0.077004
                     B     0.021358
            1  kt    A    -0.038169
                     B     0.062300
            2  kt    A    -0.133598
                     B    -0.035392
        1.0 0  kt    A    -0.077090
                     B     0.019869
            1  kt    A    -0.038920
                     B     0.062042
            2  kt    A    -0.130846
                     B    -0.034284

        plug_sites DataFrame looks like that:


                                          x  state_hist
        t   id label side plug_id
        0.0 0  kt    A    0       -0.077064          -1
                          1       -0.077064          -1
                          2       -0.077064          -1
                          3       -0.077064          -1
                     B    0        0.022936           1
                          1        0.022936           1
                          2        0.022936           1
                          3        0.022936           1
            1  kt    A    0       -0.037428          -1
                          1       -0.037428          -1
                          2       -0.037428          -1
                          3       -0.037428          -1
                     B    0        0.062572           1
                          1        0.062572           0
                          2        0.062572           1
        """
        chromosomes_dic = {}
        plugsites_dic = {}

        for i, chromosome in enumerate(chromosomes):
            ktA = chromosome.cen_A
            ktB = chromosome.cen_B

            chromosomes_dic[(i, 'kt', 'A')] = pd.DataFrame(ktA.traj, columns=['x'], index=time_index)
            chromosomes_dic[(i, 'kt', 'B')] = pd.DataFrame(ktB.traj, columns=['x'], index=time_index)

            for j, (psA, psB) in enumerate(zip(ktA.plugsites, ktB.plugsites)):
                psA_df = pd.DataFrame.from_dict({'x': psA.traj, 'state_hist': psA.state_hist})
                psA_df.index = time_index
                plugsites_dic[(i, 'kt', 'A', j)] = psA_df

                psB_df = pd.DataFrame.from_dict({'x': psB.traj, 'state_hist': psB.state_hist})
                psB_df.index = time_index
                plugsites_dic[(i, 'kt', 'B', j)] = psB_df

        kts = pd.concat(chromosomes_dic, names=['id', 'label', 'side'])
        kts = kts.reorder_levels([3, 0, 1, 2]).sort()

        plug_sites = pd.concat(plugsites_dic, names=['id', 'label', 'side', 'plug_id'])
        plug_sites = plug_sites.reorder_levels([4, 0, 1, 2, 3]).sort()
        plug_sites = plug_sites.reindex_axis(['x', 'state_hist'], axis=1)

        # Get ParamTree as Dataframe
        params = self.paramtree.to_df()
        measures = self.measuretree.to_df()

        store = pd.HDFStore(simufname)
        for dfname in df_to_save:
            store[dfname] = locals()[dfname]
        store.close()

        if verbose:
            logger.info("Simulation saved to file %s " % simufname)

    def read(self, simufname, verbose=False):
        """
        Creates a simul_spindle.Metaphase from a results.zip file.

        Parameters
        ----------
        simufname : string, optional
            The .zip file where results from the existing simulation are
            (zip which contains xml and npy)

        verbose : bool
            Set Metaphase verbose

        Returns
        -------
        :class:`~kt_simul.core.simul_spindle.Metaphase`

        """

        store = pd.HDFStore(simufname)

        param_root = SimuIO().build_tree(store['params'])
        paramtree = ParamTree(root=param_root)
        params = paramtree.relative_dic

        measure_root = SimuIO().build_tree(store['measures'])
        measuretree = ParamTree(root=measure_root,
                                adimentionalized=False)

        meta = Metaphase(paramtree=paramtree, measuretree=measuretree, verbose=False)
        KD = KinetoDynamics(params)

        spbs = store['spbs']
        KD.spbL.traj = spbs.xs('A', level='side').values.T[0]
        KD.spbR.traj = spbs.xs('B', level='side').values.T[0]

        kts = store['kts']
        plug_sites = store['plug_sites']
        store.close()

        for ch, (ch_id, ch_df) in zip(KD.chromosomes, kts.groupby(level="id")):
            ch.cen_A.traj = ch_df.xs('A', level='side')['x'].values.T
            ch.cen_B.traj = ch_df.xs('B', level='side')['x'].values.T

            for i, plugsite in enumerate(ch.cen_A.plugsites):
                ps = plug_sites.xs((ch_id, 'A', i), level=['id', 'side', 'plug_id'])
                plugsite.traj = ps['x']
                plugsite.state_hist = ps['state_hist']
                plugsite.plug_state = plugsite.state_hist[-1]

            for i, plugsite in enumerate(ch.cen_B.plugsites):
                ps = plug_sites.xs((ch_id, 'B', i), level=['id', 'side', 'plug_id'])
                plugsite.traj = ps['x']
                plugsite.state_hist = ps['state_hist']
                plugsite.plug_state = plugsite.state_hist[-1]

            ch.calc_erroneous_history()
            ch.calc_correct_history()

        meta.KD = KD
        meta.KD.simulation_done = True

        return meta

    def build_tree(self, df):
        """
        Build :class:`ElementTree` instance from :class:`pandas.DataFrame`.

        Parameters
        ----------

        df : :class:`pandas.DataFrame`
            df should come from ParamTree.to_df() method.

        Returns
        -------
        :class:`ElementTree`
        """

        root = ET.Element("parameters")
        tags = ['unit', 'description']
        for i, p in df.iterrows():
            et = ET.SubElement(root, "param")
            for key, value in p.iteritems():
                if key not in tags:
                    et.set(key, value)
                else:
                    subet = ET.SubElement(et, key)
                    subet.text = value
        return root
