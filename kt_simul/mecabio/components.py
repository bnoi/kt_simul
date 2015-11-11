import matplotlib.pylab as plt

import numpy as np
import pandas as pd

import logging
import copy

import tqdm

from . import coords, dcoords, ucoords, speed_coords
from ..utils.transformations import transformations_matrix

log = logging.getLogger(__name__)

point_cols = coords + speed_coords
coords_idxs = slice(0, 3)
speed_coords_idxs = slice(3, 6)

link_cols = dcoords + ucoords + ['length', ]
dcoords_idxs = slice(0, 3)
ucoords_idxs = slice(3, 6)
length_idx = 6


def _to_3d(df):
    df_3d = np.asarray(df).repeat(3).reshape((df.size, 3))
    return df_3d


class Structure:
    """A collection of points and links between them
    """

    def __init__(self, name):
        self.name = name

        self.points = {}
        self.point_df = np.array([])

        self.attributes = {}
        self.attributes_hist = pd.Panel()
        self._attributes_hist = []
        self.color_idx = 0

        self.point_hist = pd.Panel()
        self._point_hist = []

        self.links = {}
        self.links_idx = {}
        self.link_df = None
        self._links_idx_srce = []
        self._links_idx_trgt = []

    def add_point(self, idx=None, init_pos=None, init_speed=None, color=None):
        """Add a point to the structure
        """

        point = Point(self, idx=idx, init_pos=init_pos, init_speed=init_speed, color=color)
        return point

    def add_link(self, point_i, point_j):
        """Add a link between two point to the structure
        """

        state = np.zeros((1, len(link_cols)))

        if isinstance(self.link_df, type(None)):
            self.link_df = state
        else:
            self.link_df = np.vstack([self.link_df, state])

        link = Link(point_i, point_j, self)

        self.links[(point_i.idx, point_j.idx)] = link

        self._links_idx_srce.append(point_i.idx)
        self._links_idx_trgt.append(point_j.idx)

        self.links_idx[link.link_idx] = (point_i.idx, point_j.idx)

        return link

    @property
    def srce_idx(self):
        return np.array(self._links_idx_srce)

    @property
    def trgt_idx(self):
        return np.array(self._links_idx_trgt)

    def update_geometry(self):
        """Update link vectors according to new positions and speeds of points
        """

        self.link_df[:, dcoords_idxs] = (
            self.point_df[self.trgt_idx, coords_idxs] -
            self.point_df[self.srce_idx, coords_idxs])

        self.link_df[:, length_idx] = (self.link_df[:, dcoords_idxs]**2).sum(axis=1)**0.5

        self.link_df[:, ucoords_idxs] = (self.link_df[:, dcoords_idxs] /
                                         _to_3d(self.link_df[:, length_idx]))

        for link_idx, link in self.links.items():
            link.outer = np.outer(link.unit, link.unit)

    def register_history(self, step):
        """Save points history for re-use after the simulation
        """

        self._point_hist.append(self.point_df.copy())
        self._attributes_hist.append(copy.deepcopy(self.attributes))

    def end_history(self):
        """Various attributes conversion
        """

        self.point_hist = pd.Panel({i: pd.DataFrame(p, columns=point_cols) for i, p in enumerate(self._point_hist)})
        self._point_hist = []

        self.attributes_hist = pd.Panel({i: pd.DataFrame(p) for i, p in enumerate(self._attributes_hist)})
        self._attributes_hist = []

    def show(self):
        """Basic plot for points trajectories over time
        """

        fig = plt.figure(figsize=(14, 8))

        ax = fig.add_subplot(2, 2, 1)

        used_coords = coords.copy()
        if 'x_proj' in self.point_hist.axes[2]:
            used_coords.append('x_proj')

        for i, coord in enumerate(used_coords):
            ax = fig.add_subplot(2, 2, i+1, sharex=ax, sharey=ax)

            for n, p in self.points.items():
                if isinstance(p['color'], str):
                    ax.plot(p.traj[coord], 'o-', color=p['color'])
                else:
                    ax.plot(p.traj[coord], 'o-')

            ax.set_ylabel(coord)

        fig.tight_layout()

        return fig

    def project(self, idx_point1, idx_point2, progress=False):
        """Create a _proj column in point_hist and project all points according to an axis defined
        by two points
        """

        # Reorder Panel in DataFrame
        trajs = self.point_hist.to_frame()
        trajs = trajs.stack().unstack('minor')
        trajs = trajs.reorder_levels([1, 0]).sortlevel()
        trajs.index.names = ['time_point', 'points']

        pcoords = []
        for coord in coords:
            pcoord = '{}_proj'.format(coord)
            pcoords.append(pcoord)
            trajs[pcoord] = np.nan

        n_steps = trajs.index.get_level_values('time_point').unique().shape[0]

        iterator = tqdm.tqdm(enumerate(trajs.groupby(level='time_point')),
                             disable=not progress, total=n_steps)

        p1 = trajs.loc[pd.IndexSlice[:, idx_point1], coords].values
        p2 = trajs.loc[pd.IndexSlice[:, idx_point2], coords].values

        ref = (p1 + p2) / 2
        vec = p1 - ref

        proj_all_points = np.zeros((trajs.shape[0], len(coords)))

        for i, (t, points) in iterator:

            A = transformations_matrix(ref[i], vec[i])

            proj_points = np.dot(points[coords], A)[:, :]

            idx_start = i * len(points)
            idx_end = i * len(points) + len(points)
            proj_all_points[idx_start:idx_end] = proj_points

        trajs[pcoords] = proj_all_points

        trajs = trajs.stack().unstack('time_point')
        self.point_hist = trajs.to_panel()


class Point(object):

    def __init__(self, structure, idx=None, init_pos=None, init_speed=None, color=None):

        self.structure = structure

        if idx is None:
            self.idx = self.structure.point_df.shape[0]
        else:
            self.idx = idx
        self.idxs = slice(self.idx*3, self.idx*3 + 3)

        if init_pos is None:
            init_pos = np.zeros(len(coords))
        if init_speed is None:
            init_speed = np.zeros(len(coords))
        if isinstance(color, type(None)):
            color = np.nan

        state = np.concatenate((init_pos, init_speed)).reshape((1, len(coords)*2))

        if self.structure.point_df.shape[0] == 0:
            self.structure.point_df = state
        else:
            self.structure.point_df = np.vstack([self.structure.point_df, state])

        self['color'] = color

        self.structure.points[self.idx] = self

    @property
    def pos(self):
        return self.structure.point_df[self.idx, coords_idxs]

    @pos.setter
    def pos(self, position):
        self.structure.point_df[self.idx, coords_idxs] = position

    @property
    def speed(self):
        return self.structure.point_df[self.idx, speed_coords_idxs]

    @speed.setter
    def speed(self, speed):
        self.structure.point_df[self.idx, speed_coords_idxs] = speed

    @property
    def traj(self):
        return self.structure.point_hist.xs(self.idx, axis='major').T

    def dist(self, other):
        return np.linalg.norm(self.pos - other.pos)

    def __getitem__(self, label):
        return self.structure.attributes[label][self.idx]

    def __setitem__(self, label, value):
        if label not in self.structure.attributes.keys():
            self.structure.attributes[label] = np.zeros(self.idx + 1) * np.nan

        if isinstance(value, str):
            self.structure.attributes[label] = self.structure.attributes[label].astype('object')

        if self.idx >= self.structure.attributes[label].shape[0]:
            values = np.zeros(self.idx + 1 - self.structure.attributes[label].shape[0]) * np.nan
            self.structure.attributes[label] = np.hstack([self.structure.attributes[label], values])

        self.structure.attributes[label][self.idx] = value


class Link:

    def __init__(self, point_i, point_j, structure):

        self.structure = structure
        self.point_i = point_i
        self.point_j = point_j
        idx_i = point_i.idx
        idx_j = point_j.idx

        # Automatic link index
        self.link_idx = len(self.structure.links)

        self.idx = (idx_i, idx_j)
        self.idxs_i = point_i.idxs
        self.idxs_j = point_j.idxs

        self.idxs_ii = (self.idxs_i, self.idxs_i)
        self.idxs_ij = (self.idxs_i, self.idxs_j)
        self.idxs_jj = (self.idxs_j, self.idxs_j)
        self.idxs_ji = (self.idxs_j, self.idxs_i)

    @property
    def data(self):
        return self.structure.link_df[self.link_idx]

    @property
    def dcoords(self):
        return self.structure.link_df[self.link_idx, dcoords_idxs]

    @property
    def length(self):
        return self.structure.link_df[self.link_idx, length_idx]

    @property
    def unit(self):
        return self.structure.link_df[self.link_idx, ucoords_idxs]
