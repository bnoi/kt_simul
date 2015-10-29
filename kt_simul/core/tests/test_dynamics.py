from kt_simul.core.components import Structure, Point, Link
from kt_simul.core.dynamics import Model, dampedspring, viscous
from kt_simul.core.dynamics import spring, dashpot, linear_fv, contraction

from numpy.testing import assert_almost_equal
# This is an example test, not functional yet, should
# put some effort in implementing for all forces

def test_spring():

    sprg = Structure('simple spring')
    p0 = Point(0, sprg)
    sprg.add_point(p0, pos0=[0, 0, 0])
    p1 = Point(1, sprg)
    sprg.add_point(p1, pos0=[2, 0, 0])
    lnk = sprg.add_link(p0, p1)
    sprg.update_geometry()
    m = Model(sprg)
    viscous(m, p0, 1)
    viscous(m, p1, 1)
    sprg.point_hist = pd.Panel({0: sprg.point_df})

    def model_update(step):

        m.Bvect *= 0
        spring(m, lnk, 0.1, 1)
        sprg.register_history(step)

    for i in range(30):
        m.solve()
        model_update(i)

    assert_almost_equal(p0.dist(p1), 1, decimal=2)
