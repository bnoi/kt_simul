from kt_simul.mecabio import Structure
from kt_simul.mecabio import Point

from numpy.testing import assert_almost_equal


def test_structure():
    """Test points, links and attributes data structure
    """

    struct = Structure('')
    p0 = struct.add_point(0, init_pos=[2, 0, 0], color="black")

    assert_almost_equal(struct.point_df[0], [2., 0., 0., 0., 0., 0.])
    assert struct.attributes_df.loc[0, 'color'] == 'black'

    struct = Structure('')
    p0 = Point(struct, 0, init_pos=[2, 0, 0], color="black")
    p1 = Point(struct, init_pos=[-2, 0, 0], color="red")

    struct.add_link(p0, p1)

    assert struct.point_df.shape == (2, 6)
    assert struct.link_df.shape == (1, 7)
