from kt_simul.core import Metaphase
from kt_simul.core import parameters


def test_metaphase():
    """Only check wether basic spindle simulation end without exception
    """

    params = parameters.get_default_params()
    params.loc['dt', 'value'] = 1
    params.loc['span', 'value'] = 10
    params.loc['t_A', 'value'] = 10
    measures = parameters.get_default_measures()

    meta = Metaphase(verbose=False,
                     params=params,
                     measures=measures,
                     initial_plug='null',
                     keep_same_random_seed=False,
                     force_parameters=[])

    meta.simul(progress=False)
    meta.project(progress=False)


def test_dt_coherence():
    """
    """
    params = parameters.get_default_params()
    params.loc['dt', 'value'] = 0.1
    params.loc['span', 'value'] = 2
    params.loc['t_A', 'value'] = 2

    meta = Metaphase(verbose=False, params=params)
    meta.simul()

    excepted_n_steps = (params.loc['span', 'value'] / params.loc['dt', 'value'])

    assert meta.spindle.point_hist.shape[0] == excepted_n_steps
