import numpy as np

from . import coords


class Model:

    def __init__(self, structure, dt=1):

        self.structure = structure
        self.dt = dt
        self.n_points = len(self.structure.points)
        self.Amat = np.zeros((self.n_points * 3, self.n_points * 3))
        self.Bvect = np.zeros(self.n_points * 3)

        for point in self.structure.points.values():
            viscous(self, point, 1)

    def solve(self):
        speeds = np.linalg.solve(self.Amat, -self.Bvect)
        for i, c in enumerate(coords):
            self.structure.point_df['v'+c] = speeds[i::3]
            self.structure.point_df[c] += speeds[i::3] * self.dt
        self.structure.update_geometry()


def viscous(model, point, mu):
    ''' Viscous drag with coeff mu on point '''
    model.Amat[point.idxs, point.idxs] = -np.eye(3) * mu


def dashpot(model, link, mu):

    block = mu * link.outer
    model.Amat[link.idxs_ii] -= block
    model.Amat[link.idxs_ij] += block
    model.Amat[link.idxs_jj] -= block
    model.Amat[link.idxs_ji] += block


def spring(model, link, kappa, d_eq=0):

    idx_i, idx_j = link.idx

    F = kappa * (link.length - d_eq) * link.unit
    model.Bvect[link.idxs_i] += F
    model.Bvect[link.idxs_j] -= F


def dampedspring(model, link, mu, kappa, d_eq=0):

    dashpot(model, link, mu)
    spring(model, link, kappa, d_eq)


def contraction(model, link, F):

    model.Bvect[link.idxs_i] += F * link.unit
    model.Bvect[link.idxs_j] -= F * link.unit


def linear_fv(model, link, F_stall, v_max, gamma=-1):

    dashpot(model, link, F_stall / v_max)
    contraction(model, link, -gamma * F_stall)
