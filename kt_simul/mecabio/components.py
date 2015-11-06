import numpy as np
import pandas as pd
import logging

from . import coords, dcoords, ucoords, speed_coords

log = logging.getLogger(__name__)

point_cols = coords + speed_coords
link_cols = dcoords + ucoords + ['length', ]


def _to_3d(df):
    df_3d = np.asarray(df).repeat(3).reshape((df.size, 3))
    return df_3d


class Structure:
    '''
    A collection of points and links between them
    '''
    def __init__(self, name):
        self.name = name
        self.points = {}
        self.links = {}
        self.point_df = pd.DataFrame(columns=point_cols)
        _link_idx = pd.MultiIndex(levels=[[], []],
                                  labels=[[], []],
                                  names=['srce', 'trgt'])
        self.link_df = pd.DataFrame(columns=link_cols, index=_link_idx)
        self.point_hist = pd.Panel()

    def add_point(self, point, pos0=None, speed0=None):

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

        self.link_df[dcoords] = (
            self.point_df[coords].loc[self.trgt_idx].values -
            self.point_df[coords].loc[self.srce_idx].values)

        self.link_df['length'] = (self.link_df[dcoords]**2).sum(axis=1)**0.5

        self.link_df[ucoords] = (self.link_df[dcoords] /
                                 _to_3d(self.link_df['length']))

        for idxs, link_values in self.link_df[ucoords].iterrows():
            link = self.links[idxs]
            link.outer = np.outer(link_values, link_values)

        self.link_df.sortlevel(inplace=True)
        self.point_df.sortlevel(inplace=True)

    def register_history(self, step):
        self.point_hist[step] = self.point_df


class Point:

    def __init__(self, idx, structure):

        self.idx = idx
        self.idxs = slice(idx*3, idx*3 + 3)
        self.structure = structure

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
