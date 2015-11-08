from kt_simul.core import Metaphase
from kt_simul.core import load_metaphase
from kt_simul.core import parameters


def test_metaphase():
    """Only check is basic spindle simulation end without exception
    """

    params = parameters.get_default_params()
    params.loc['dt', 'value'] = 1
    params.loc['span', 'value'] = 10
    params.loc['t_A', 'value'] = 10

    params.loc['d_alpha', 'value'] = 0.05  # 0.05
    params.loc['k_a', 'value'] = 0.1  # 0.06
    params.loc['orientation', 'value'] = 1  # 1

    params.loc['N', 'value'] = 1  # 3
    params.loc['Mk', 'value'] = 1  # 3
    params.loc['L0', 'value'] = 0.3  # 0.3

    params.loc['ldep', 'value'] = 0  # 0.2

    measures = parameters.get_default_measures()

    meta = Metaphase(verbose=False,
                     params=params,
                     measures=measures,
                     initial_plug='null',
                     keep_same_random_seed=False,
                     force_parameters=[])

    meta.simul(progress=False)
    meta.project(progress=False)
