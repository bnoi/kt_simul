from kt_simul.mecabio import Structure
from kt_simul.mecabio import Point

from kt_simul.mecabio import Model

from kt_simul.mecabio import spring
from kt_simul.mecabio import dampedspring
from kt_simul.mecabio import viscous
from kt_simul.mecabio import dashpot
from kt_simul.mecabio import linear_fv
from kt_simul.mecabio import contraction

from numpy.testing import assert_almost_equal


def test_spring():

    sprg = Structure('')
    p0 = Point(0, sprg)
    sprg.add_point(p0, pos0=[0, 0, 0])
    p1 = Point(1, sprg)
    sprg.add_point(p1, pos0=[2, 0, 0])
    lnk = sprg.add_link(p0, p1)
    sprg.update_geometry()
    m = Model(sprg)

    def model_update(step):

        m.Bvect *= 0
        spring(m, lnk, 0.1, 1)
        sprg.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 1, decimal=2)


def test_dampedspring():

    sprg = Structure('')
    p0 = Point(0, sprg)
    sprg.add_point(p0, pos0=[2, 0, 0])
    p1 = Point(1, sprg)
    sprg.add_point(p1, pos0=[-2, 0, 0])
    lnk = sprg.add_link(p0, p1)
    sprg.update_geometry()
    m = Model(sprg)

    def model_update(step):

        m.Bvect *= 0
        dampedspring(m, lnk, 0.1, 0.1, 1)
        sprg.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 1.44, decimal=2)


def test_viscous():

    sprg = Structure('')
    p0 = Point(0, sprg)
    sprg.add_point(p0, pos0=[2, 0, 0])
    p1 = Point(1, sprg)
    sprg.add_point(p1, pos0=[-2, 0, 0])
    lnk = sprg.add_link(p0, p1)
    sprg.update_geometry()
    m = Model(sprg)

    def model_update(step):

        m.Bvect *= 0
        contraction(m, lnk, 1)
        viscous(m, p0, 100)
        viscous(m, p1, 100)
        sprg.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 3.42, decimal=2)


def test_dashpot():

    sprg = Structure('')
    p0 = Point(0, sprg)
    sprg.add_point(p0, pos0=[2, 0, 0])
    p1 = Point(1, sprg)
    sprg.add_point(p1, pos0=[-2, 0, 0])
    lnk = sprg.add_link(p0, p1)
    sprg.update_geometry()
    m = Model(sprg)

    def model_update(step):

        m.Bvect *= 0
        dashpot(m, lnk, 1)
        contraction(m, lnk, 1)
        sprg.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 0.64, decimal=2)


def test_linear_fv():
    sprg = Structure('')
    p0 = Point(0, sprg)
    sprg.add_point(p0, pos0=[2, 0, 0])
    p1 = Point(1, sprg)
    sprg.add_point(p1, pos0=[-2, 0, 0])
    lnk = sprg.add_link(p0, p1)
    sprg.update_geometry()
    m = Model(sprg)

    def model_update(step):

        m.Bvect *= 0
        linear_fv(m, lnk, 1, 1, 1)
        sprg.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 7.36, decimal=2)


def test_contraction():
    sprg = Structure('')
    p0 = Point(0, sprg)
    sprg.add_point(p0, pos0=[2, 0, 0])
    p1 = Point(1, sprg)
    sprg.add_point(p1, pos0=[-2, 0, 0])
    lnk = sprg.add_link(p0, p1)
    sprg.update_geometry()
    m = Model(sprg)

    def model_update(step):

        m.Bvect *= 0
        contraction(m, lnk, 0.05)
        sprg.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 1.10, decimal=2)
