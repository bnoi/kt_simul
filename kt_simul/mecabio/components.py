import matplotlib.pylab as plt

import numpy as np
import pandas as pd

import logging

import tqdm

from . import coords, dcoords, ucoords, speed_coords
from ..utils.transformations import transformations_matrix

log = logging.getLogger(__name__)

point_cols = coords + speed_coords
link_cols = dcoords + ucoords + ['length', ]


def _to_3d(df):
    df_3d = np.asarray(df).repeat(3).reshape((df.size, 3))
    return df_3d


class Structure:
    """A collection of points and links between them
    """

    def __init__(self, name):
        self.name = name
        self.points = {}
        self.links = {}
        self.point_df = pd.DataFrame(columns=point_cols)
        _link_idx = pd.MultiIndex(levels=[[], []],
                                  labels=[[], []],
                                  names=['srce', 'trgt'])
        self.link_df = pd.DataFrame(columns=link_cols, index=_link_idx)
        self.point_hist = None

    def add_point(self, point, pos0=None, speed0=None):
        """Add a point to the structure
        """

        if point.idx is None:
            point.idx = self.point_df.shape[0]
        if pos0 is None:
            pos0 = np.zeros(len(coords))
        if speed0 is None:
            speed0 = np.zeros(len(coords))

        state = np.concatenate((pos0, speed0)).reshape((1, len(coords)*2))
        self.point_df = self.point_df.append(pd.DataFrame(index=[point.idx, ],
                                                          data=state,
                                                          columns=point_cols))
        self.points[point.idx] = point

    def add_link(self, point_i, point_j):
        """Add a link between two point to the structure
        """

        state = np.zeros((1, self.link_df.shape[1]))
        self.link_df = self.link_df.append(
            pd.DataFrame(index=[(point_i.idx, point_j.idx)],
                         data=state,
                         columns=self.link_df.columns))
        link = Link(point_i, point_j, self)
        self.links[(point_i.idx, point_j.idx)] = link
        return link

    @property
    def srce_idx(self):
        return self.link_df.index.get_level_values('srce')

    @property
    def trgt_idx(self):
        return self.link_df.index.get_level_values('trgt')

    def update_geometry(self):
        """Update link vectors according to new positions and speeds of points
        """

        self.link_df[dcoords] = (
            self.point_df[coords].loc[self.trgt_idx].values -
            self.point_df[coords].loc[self.srce_idx].values)

        self.link_df['length'] = (self.link_df[dcoords]**2).sum(axis=1)**0.5

        self.link_df[ucoords] = (self.link_df[dcoords] /
                                 _to_3d(self.link_df['length']))

        for idxs, link_values in self.link_df[ucoords].iterrows():
            link = self.links[idxs]
            link.outer = np.outer(link_values, link_values)

    def register_history(self, step):
        """Save points history for re-use after the simulation
        """

        if isinstance(self.point_hist, type(None)):
            self.point_hist = pd.Panel({step: self.point_df})
        else:
            self.point_hist[step] = self.point_df

    def show(self):
        """Basic plot for points trajectories over time
        """

        fig = plt.figure(figsize=(18, 12))

        ax = fig.add_subplot(len(coords)//2+1, 2, 1)

        used_coords = coords.copy()
        if 'x_proj' in self.point_hist.axes[2]:
            used_coords.append('x_proj')

        for i, coord in enumerate(used_coords):
            ax = fig.add_subplot(len(used_coords)//2+1, 2, i+1, sharex=ax, sharey=ax)

            for n, p in self.points.items():
                ax.plot(p.traj[coord], 'o-', color=p.color)

            ax.set_ylabel(coord)

        return fig

    def project(self, idx_point1, idx_point2, progress=False):
        """Create a _proj column in point_hist and project all points according to an axis defined
        by two points
        """

        # Reorder Panel in DataFrame
        trajs = self.point_hist.to_frame()
        trajs = trajs.stack().unstack('minor')
        trajs = trajs.reorder_levels([1, 0]).sortlevel()
        trajs.index.names = ['times', 'points']

        for coord in coords:
            trajs['{}_proj'.format(coord)] = np.nan

        duration = trajs.index.get_level_values('times').max()

        iterator = tqdm.tqdm(enumerate(trajs.groupby(level='times')),
                             disable=not progress, total=duration)

        for i, (t, points) in iterator:

            p1 = points.loc[t, idx_point1][coords]
            p2 = points.loc[t, idx_point2][coords]

            ref = (p1 + p2) / 2
            vec = p1 - ref

            A = transformations_matrix(ref, vec)

            projected_points = np.dot(points[coords], A)[:, :]

            for i, coord in enumerate(coords):
                trajs.loc[t, '{}_proj'.format(coord)] = projected_points[:, i]

        trajs = trajs.stack().unstack('times')
        self.point_hist = trajs.to_panel()


class Point:

    def __init__(self, idx, structure, color=None):

        self.idx = idx
        self.idxs = slice(idx*3, idx*3 + 3)
        self.structure = structure
        self._color = color

    @property
    def pos(self):
        return self.structure.point_df[coords].values[self.idx]

    @pos.setter
    def pos(self, position):
        self.structure.point_df.loc[self.idx, coords] = position

    @property
    def speed(self):
        return self.structure.point_df.loc[self.idx, speed_coords]

    @speed.setter
    def speed(self, speed):
        self.structure.point_df.loc[self.idx, speed_coords] = speed

    @property
    def traj(self):
        return self.structure.point_hist.xs(self.idx, axis='major').T

    def dist(self, other):
        return np.linalg.norm(self.pos - other.pos)

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, color):
        self._color = color


class Link:

    def __init__(self, point_i, point_j, structure):

        self.structure = structure
        self.point_i = point_i
        self.point_j = point_j
        idx_i = point_i.idx
        idx_j = point_j.idx

        self.idx = (idx_i, idx_j)
        self.idxs_i = point_i.idxs
        self.idxs_j = point_j.idxs

        self.idxs_ii = (self.idxs_i, self.idxs_i)
        self.idxs_ij = (self.idxs_i, self.idxs_j)
        self.idxs_jj = (self.idxs_j, self.idxs_j)
        self.idxs_ji = (self.idxs_j, self.idxs_i)

    @property
    def data(self):
        return self.structure.link_df.loc[self.idx]

    @property
    def dcoords(self):
        return self.structure.link_df.loc[self.idx, dcoords]

    @property
    def length(self):
        return self.structure.link_df.loc[self.idx, 'length']

    @property
    def speed(self):
        return self.structure.link_df.loc[self.idx, dspeeds]

    @property
    def unit(self):
        return self.structure.link_df.loc[self.idx][ucoords]
