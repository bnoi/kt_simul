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


def test_datastructure():
    """Test points, links and attributes data structure
    """

    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[2, 0, 0], color="black")

    assert_almost_equal(struct.point_df.values[0], [2., 0., 0., 0., 0., 0.])
    assert struct.attributes_df.loc[0, 'color'] == 'black'

    struct = Structure('')
    p0 = Point(struct, 0, init_pos=[2, 0, 0], color="black")
    p1 = Point(struct, init_pos=[-2, 0, 0], color="red")

    struct.add_link(p0, p1)

    assert struct.point_df.index[0] == 0
    assert struct.link_df.index[0] == (0, 1)
    assert struct.link_df.shape == (1, 7)


def test_spring():

    sprg = Structure('')
    p0 = sprg.add_point(0, init_pos=[0, 0, 0])
    p1 = sprg.add_point(1, init_pos=[2, 0, 0])
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
    p0 = sprg.add_point(0, init_pos=[0, 0, 0])
    p1 = sprg.add_point(1, init_pos=[2, 0, 0])
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

    assert_almost_equal(p0.dist(p1), 1.15, decimal=2)


def test_viscous():

    sprg = Structure('')
    p0 = sprg.add_point(0, init_pos=[0, 0, 0])
    p1 = sprg.add_point(1, init_pos=[2, 0, 0])
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

    assert_almost_equal(p0.dist(p1), 1.42, decimal=2)


def test_dashpot():

    sprg = Structure('')
    p0 = sprg.add_point(0, init_pos=[0, 0, 0])
    p1 = sprg.add_point(1, init_pos=[2, 0, 0])
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

    assert_almost_equal(p0.dist(p1), 0.00, decimal=2)


def test_linear_fv():
    sprg = Structure('')
    p0 = sprg.add_point(0, init_pos=[0, 0, 0])
    p1 = sprg.add_point(1, init_pos=[2, 0, 0])
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

    assert_almost_equal(p0.dist(p1), 5.36, decimal=2)


def test_contraction():
    sprg = Structure('')
    p0 = sprg.add_point(0, init_pos=[0, 0, 0])
    p1 = sprg.add_point(1, init_pos=[2, 0, 0])
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

    assert_almost_equal(p0.dist(p1), 0.10, decimal=2)
