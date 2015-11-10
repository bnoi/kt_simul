import pandas as pd

from kt_simul.mecabio import Structure

from kt_simul.mecabio import Model

from kt_simul.mecabio import spring
from kt_simul.mecabio import dampedspring
from kt_simul.mecabio import viscous
from kt_simul.mecabio import dashpot
from kt_simul.mecabio import linear_fv
from kt_simul.mecabio import contraction
from kt_simul.mecabio.dynamics import PhysicsException

from numpy.testing import assert_almost_equal

from nose.tools import raises


def test_history():
    """
    """
    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[0, 0, 0])
    p1 = struct.add_point(1, init_pos=[2, 0, 0])
    lnk = struct.add_link(p0, p1)
    struct.update_geometry()
    m = Model(struct)

    viscous(m, p0, 1)
    viscous(m, p1, 1)

    def model_update(step):

        m.Bvect *= 0
        spring(m, lnk, 0.1, 1)
        struct.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    struct.end_history()

    assert struct.point_hist.shape == (30, 2, 6)
    assert isinstance(struct.point_hist, pd.Panel)


def test_projection():
    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[0, 0, 0])
    p1 = struct.add_point(1, init_pos=[2, 0, 0])
    lnk = struct.add_link(p0, p1)
    struct.update_geometry()
    m = Model(struct)

    viscous(m, p0, 1)
    viscous(m, p1, 1)

    def model_update(step):

        m.Bvect *= 0
        spring(m, lnk, 0.1, 1)
        struct.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    struct.end_history()
    struct.project(0, 1)

    assert set(['x_proj', 'y_proj', 'z_proj']).issubset(set(struct.point_hist.loc[0].columns))
    assert_almost_equal(struct.point_hist.loc[0]['x_proj'].values, [0, -2])


def test_spring():

    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[0, 0, 0])
    p1 = struct.add_point(1, init_pos=[2, 0, 0])
    lnk = struct.add_link(p0, p1)
    struct.update_geometry()
    m = Model(struct)

    viscous(m, p0, 1)
    viscous(m, p1, 1)

    def model_update(step):

        m.Bvect *= 0
        spring(m, lnk, 0.1, 1)
        struct.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 1, decimal=2)


def test_dampedspring():

    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[0, 0, 0])
    p1 = struct.add_point(1, init_pos=[2, 0, 0])
    lnk = struct.add_link(p0, p1)
    struct.update_geometry()
    m = Model(struct)

    viscous(m, p0, 1)
    viscous(m, p1, 1)

    def model_update(step):

        m.Bvect *= 0
        dampedspring(m, lnk, 0.1, 0.1, 1)
        struct.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 1.15, decimal=2)


def test_viscous():

    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[0, 0, 0])
    p1 = struct.add_point(1, init_pos=[2, 0, 0])
    lnk = struct.add_link(p0, p1)
    struct.update_geometry()
    m = Model(struct)

    viscous(m, p0, 1)
    viscous(m, p1, 1)

    def model_update(step):

        m.Bvect *= 0
        contraction(m, lnk, 1)
        viscous(m, p0, 100)
        viscous(m, p1, 100)
        struct.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 1.42, decimal=2)


def test_dashpot():

    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[0, 0, 0])
    p1 = struct.add_point(1, init_pos=[2, 0, 0])
    lnk = struct.add_link(p0, p1)
    struct.update_geometry()
    m = Model(struct)

    viscous(m, p0, 1)
    viscous(m, p1, 1)

    def model_update(step):

        m.Bvect *= 0
        dashpot(m, lnk, 1)
        contraction(m, lnk, 1)
        struct.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 0.00, decimal=2)


def test_linear_fv():
    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[0, 0, 0])
    p1 = struct.add_point(1, init_pos=[2, 0, 0])
    lnk = struct.add_link(p0, p1)
    struct.update_geometry()
    m = Model(struct)

    viscous(m, p0, 1)
    viscous(m, p1, 1)

    def model_update(step):

        m.Bvect *= 0
        linear_fv(m, lnk, 1, 1, 1)
        struct.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 5.36, decimal=2)


def test_contraction():
    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[0, 0, 0])
    p1 = struct.add_point(1, init_pos=[2, 0, 0])
    lnk = struct.add_link(p0, p1)
    struct.update_geometry()
    m = Model(struct)

    viscous(m, p0, 1)
    viscous(m, p1, 1)

    def model_update(step):

        m.Bvect *= 0
        contraction(m, lnk, 0.05)
        struct.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 0.10, decimal=2)


@raises(PhysicsException)
def test_viscous_points():
    """
    """
    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[0, 0, 0])
    p1 = struct.add_point(1, init_pos=[2, 0, 0])
    lnk = struct.add_link(p0, p1)
    struct.update_geometry()
    m = Model(struct)

    viscous(m, p0, 1)

    def model_update(step):

        m.Bvect *= 0
        spring(m, lnk, 0.1, 1)
        struct.register_history(step)

    for i in range(5):
        m.solve()
        model_update(i)
