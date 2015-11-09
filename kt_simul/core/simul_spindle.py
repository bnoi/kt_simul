"""
This module provides the core simulation functionalities.

See Gay et al. J. Cell Biol., 2012 http://dx.doi.org/10.1083/jcb.201107124
The original framework was adapted from:
Civelekoglu-Scholey et al. Biophys. J. 90(11), 2006
http://dx.doi.org/10.1529/biophysj.105.078691
"""

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import logging
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

import tqdm

from ..core.components import Spindle
from ..core.dynamics import SpindleModel
from ..core import parameters

log = logging.getLogger(__name__)

__all__ = ["Metaphase"]


class Metaphase(object):

    """
    An instance of the Metaphase class is a wrapper around
    the whole simulation.

    Parameters
    ----------

    params : dict
        The params contains the parameters for the simulation
        if params is None, the parameters are read
        from the file paramfile. Defaults to None.

    measures : dict
        The measures contains the observed characteristics
        of the mitosis e.g. metaphase spindle elongation rate, etc.
        if measures is None, the measures are read from the file
        indicated by the measurefile argument. Defaults to None.

    initial_plug : string or None
        Defines globally the initial attachment states. This argument can have
        the following values:
            - 'null': all kinetochores are detached
            - 'amphitelic': all chromosmes are amphitelic
            - 'random': all attachement site can be bound to either pole or
            detected with equal prob.
            - 'monotelic': right kinetochores are attached to the same pole,
            left ones are detached
            - 'syntelic' : all kinetochores are attached to the same pole

    keep_same_random_seed : bool
        To launch simulations with same random state

    """

    RANDOM_STATE = None

    def __init__(self,  params=None, measures=None,
                 initial_plug='random', verbose=False,
                 keep_same_random_seed=False,
                 force_parameters=[], name='sim'):

        # Enable or disable log console
        self.verbose = verbose

        if not self.verbose:
            log.disabled = True
        else:
            log.disabled = False

        if params is None:
            self.params = parameters.get_default_params()
        else:
            self.params = params

        if measures is None:
            self.measures = parameters.get_default_measures()
        else:
            self.measures = measures

        self.original_params = self.params.copy()
        self.original_measures = self.measures.copy()

        self.measures = self.measures['value'].to_dict()
        self.params = self.params['value'].to_dict()

        # Parameters reduction
        self.params = parameters.reduce_params(self.params, self.measures,
                                               force_parameters=force_parameters)
        self.params = parameters.adimentionalize(self.params, self.original_params)

        log.info('Parameters loaded')

        if keep_same_random_seed:
            self.prng = self.__class__.get_random_state()
        else:
            self.prng = np.random.RandomState()

        self.spindle = Spindle(name, self.params,
                               initial_plug=initial_plug,
                               prng=self.prng)

        self.model = SpindleModel(self.spindle)

        dt = self.params['dt']
        duration = self.params['span']

        self.num_steps = int(duration / dt)
        self.model.anaphase = False
        self.timelapse = np.arange(0, duration, dt)
        self.report = []
        self.delay = -1
        self.observations = {}
        self.analysis = {}

        log.info('Simulation initialized')

    @classmethod
    def get_random_state(cls):
        """
        """
        prng = np.random.RandomState()
        if not cls.RANDOM_STATE:
            cls.RANDOM_STATE = prng.get_state()
        else:
            prng.set_state(cls.RANDOM_STATE)
        return prng

    @property
    def time(self):
        return np.arange(0, self.params['span'], self.model.dt)

    @property
    def time_anaphase(self):
        return self.analysis['real_t_A']

    @property
    def index_anaphase(self):
        return np.argwhere(self.time == self.time_anaphase)[0][0]

    def index_anaphase_before(self, t_shift):
        return np.argwhere(self.time == (self.time_anaphase - t_shift))[0][0]

    def simul(self, ablat=None, ablat_pos=0., progress=False):
        """
        The simulation main loop.

        Parameters
        ----------

        ablat : float, optional
            Timepoint at which ablation takes place. If None (default)
            no ablation is performed.

        """

        kappa_c = self.model.params['kappa_c']

        log.info('Running simulation')

        self.analysis['real_t_A'] = self.model.params['t_A']

        for time_point in tqdm.tqdm(range(0, self.num_steps),
                                    total=self.num_steps,
                                    disable=not progress):

            # Ablation test
            if ablat == time_point:
                log.info("Performing ablation")
                self._ablation(time_point, pos=ablat_pos)

            # Anaphase transition ?
            if self._anaphase_test(time_point):
                log.info("Anaphase onset at {} / {}".format(time_point, self.num_steps))

            self.model.one_step(time_point)

        log.info('Simulation done')

        self.model.params['kappa_c'] = kappa_c

        for ch in self.spindle.chromosomes:
            ch.calc_correct_history()
            ch.calc_erroneous_history()
            ch.cen_A.calc_toa()
            ch.cen_B.calc_toa()

    def save(self, fname, save_tree=True):
        """
        """

        save_params = True

        with pd.HDFStore(fname) as store:
            store["point_hist"] = self.spindle.point_hist
            store["link_df"] = pd.DataFrame(self.spindle.link_df)
            store["attributes_df"] = self.spindle.attributes_df

            store["general_params"] = pd.Series(dict(initial_plug=self.spindle.initial_plug))
            store["analysis"] = pd.Series(self.analysis)

            if save_params:
                store['params'] = pd.Series(self.params)
                store['measures'] = pd.Series(self.measures)
                store['original_params'] = self.original_params
                store['original_measures'] = self.original_measures

        log.info("Simu saved to {}".format(fname))

    def project(self, progress=False):
        """Project all points coordinates on pole-pole axes
        """

        self.spindle.project(self.spindle.spbL.idx,
                             self.spindle.spbR.idx,
                             progress=progress)

    def _anaphase_test(self, time_point):
        """
        At anaphase onset, set the cohesin spring constant to 0 and
        self.model.anaphase to True.

        Returns
        -------
        bool
            True if anaphase has been executed.

        """

        t_A = int(self.model.params['t_A'])
        dt = self.model.params['dt']
        t = time_point * dt

        if self.model.anaphase:
            return True

        if t >= t_A and self._plug_checkpoint():
            self.delay = t - t_A
            self.analysis['real_t_A'] = t

            # Then we just get rid of cohesin
            self.model.params['kappa_c'] = 0
            nb_mero = self._mero_checkpoint()
            if nb_mero:
                s = ("There were %d merotelic MT at anaphase onset"
                     % nb_mero)
                self.report.append(s)
            self.model.anaphase = True
            return True
        return False

    def _ablation(self, time_point, pos=None, radius=0.2):
        """
        Simulates a laser ablation: detaches all the kinetochores
        and sets the midzone stall force Fmz to 0

        Parameters
        ----------

        pos : float or None, optional
            Position of the laser beam within the spindle

        """
        if pos is None:
            pos = self.spindle.spbR.pos
        if not self.spindle.spbL.pos.x <= pos[0] <= self.spindle.spbR.pos.x:
            log.warning('Missed shot, same player play again!')
            return
        self.model.params['Fmz'] = 0.
        self.model.params['k_a'] = 0.
        self.model.params['k_d0'] = 0.
        self.model.A0_mat = self.model.time_invariantA()

        for plugsite in self.spindle.all_plugsites():
            if pos < plugsite.pos.x and plugsite.plug_state == - 1:
                plugsite.plug_state = 0
            elif pos > plugsite.pos and plugsite.plug_state == 1:
                plugsite.plug_state = 0

    def _plug_checkpoint(self):
        """If the spindle assembly checkpoint is active, returns True
        if all chromosomes are plugged by at least one kMT, False
        otherwise.
        """
        sac = self.model.params['sac']
        if sac == 0:
            return True

        for ch in self.spindle.chromosomes:
            if not ch.cen_A.is_attached() or not ch.cen_B.is_attached():
                return False
        return True

    def _mero_checkpoint(self):
        """
        Returns
        -------

        int
            The total number of merotellic kT
        """
        nb_mero = 0
        for ch in self.spindle.chromosomes:
            if np.any(ch.erroneous()):
                nb_mero += sum(ch.erroneous())
        return nb_mero

    def _mplate_checkpoint(self):
        """
        Returns
        -------
        bool
            Returns True if each kinetochore is in the proper half
            of the spindle
        """
        for ch in self.spindle.chromosomes:
            ktR = ch.cen_A.pos.x
            ktL = ch.cen_B.pos.x
            if min(ktR, ktL) <= 0 and max(ktR, ktL) >= 0:
                return True
        return True

    def show(self):
        """Quickly show kymograph of current simulation.
        """

        fig = plt.figure(figsize=(18, 12))

        coords = ['x', 'y', 'z']
        if 'x_proj' in self.spindle.point_hist.axes[2]:
            coords.append('x_proj')

        max_pos = max(self.spindle.spbL.traj[coords].values.flatten().max(),
                      self.spindle.spbR.traj[coords].values.flatten().max())

        ax = fig.add_subplot(2, 2, 1)
        for i, coord in enumerate(coords):

            ax = fig.add_subplot(2, 2, i+1, sharex=ax, sharey=ax)

            ax.plot(self.time, self.spindle.spbL.traj[coord],
                    color=self.spindle.spbL['color'], lw=2)
            ax.plot(self.time, self.spindle.spbR.traj[coord],
                    color=self.spindle.spbL['color'], lw=2)

            ax.set_xlabel("Time (s)", fontsize=18)
            ax.set_ylabel("{} (um)".format(coord), fontsize=18)

            for ch in self.spindle.chromosomes:
                ax.plot(self.time, ch.cen_A.traj[coord], color=ch.cen_A['color'], lw=2)
                ax.plot(self.time, ch.cen_B.traj[coord], color=ch.cen_B['color'], lw=2)

            ax.set_ylim(-max_pos, max_pos)
            ax.set_xlim(0, self.time[-1])

            ax.axvline(self.time_anaphase, lw=1, alpha=0.8, color="black")

        fig.tight_layout()

        return fig

    def get_attachment_vector(self):
        """Get attachment states for all chromosomes and all timepoints.
        """

        att = []

        for ch in self.spindle.chromosomes:
            ch.calc_correct_history()
            ch.calc_erroneous_history()

            # Return attachment history for all timepoints
            c_hist = ch.correct_history
            e_hist = ch.erroneous_history
            state = [self.get_attachment(np.concatenate((c, e)))
                     for c, e in zip(c_hist, e_hist)]

            att.append(state)

        return np.array(att).T

    def get_attachment(self, state):
        """
        Get attachment state name according to biological classification:
            - amphitelic returns 0
            - monotelic returns 1
            - syntelic returns 2
            - merotelic returns 3
            - unattached returns 4
            - error means no attachment were found (issue with code) returns 5
        """

        def amphitelic(state):
            return (state[0] > 0 and state[1] > 0 and
                    state[2] == 0 and state[3] == 0)

        def monotelic(state):
            return (state[0] > 0 and state[1] == 0 and
                    state[2] == 0 and state[3] == 0)

        def syntelic(state):
            return (state[0] > 0 and state[1] == 0 and
                    state[2] == 0 and state[3] > 0)

        def merotelic(state):
            return (state[0] > 0 and state[2] > 0)

        def unattached(state):
            return (state[0] == 0 and state[1] == 0 and
                    state[2] == 0 and state[3] == 0)

        permuted_state = [state[1], state[0], state[3], state[2]]

        if amphitelic(state) or amphitelic(permuted_state):
            return 0
        elif monotelic(state) or monotelic(permuted_state):
            return 1
        elif syntelic(state) or syntelic(permuted_state):
            return 2
        elif merotelic(state) or merotelic(permuted_state):
            return 3
        elif unattached(state) or unattached(permuted_state):
            return 4
        else:
            return 5

    def get_attachment_names(self):
        """
        """

        names = ['amphitelic', 'monotelic', 'syntelic',
                 'merotelic', 'unattached', 'error']
        return names

    def __del__(self):
        """
        """
        if hasattr(self, 'spindle'):
            del self.spindle
        if hasattr(self, 'model'):
            del self.model
        if hasattr(self, 'analysis'):
            del self.analysis
        if hasattr(self, 'measures'):
            del self.measures
        if hasattr(self, 'params'):
            del self.params


def load_metaphase(fname, params=None, measures=None, verbose=True):
    """
    """

    no_params = False
    if not params or not measures:
        no_params = True

    with pd.HDFStore(fname) as store:
        if ('params' not in store) or ('measures' not in store) and no_params:
            raise Exception("Can't load simulation without params and measures")

        if no_params:
            params = store['original_params']
            measures = store['original_measures']

        point_hist = store["point_hist"]
        link_df = store["link_df"].values
        attributes_df = store["attributes_df"]

        general_params = store["general_params"]
        analysis = store["analysis"]

    meta = Metaphase(verbose=verbose,
                     params=params,
                     measures=measures,
                     initial_plug=general_params["initial_plug"])

    meta.spindle.link_df = link_df
    meta.spindle.point_hist = point_hist
    meta.spindle.attributes_df = attributes_df

    meta.analysis = analysis

    return meta
