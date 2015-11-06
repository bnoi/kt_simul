"""
This packages contains the kinetochore dynamics simulation
"""

import logging

from .utils.color_system import color

__version__ = "2.0"


def in_ipython():
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True

if in_ipython():
    logformat = '%(asctime)s' + ':'
    logformat += '%(levelname)s' + ':'
    # logformat += '%(name)s' + ':'
    # logformat += '%(funcName)s' + ': '
    logformat += ' %(message)s'
else:
    logformat = color('%(asctime)s', 'BLUE') + ':'
    logformat += color('%(levelname)s', 'RED') + ':'
    # logformat += color('%(name)s', 'YELLOW') + ':'
    # logformat += color('%(funcName)s', 'GREEN') + ': '
    logformat += color(' %(message)s', 'ENDC')

log = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(logformat, "%Y-%m-%d %H:%M:%S")
handler.setFormatter(formatter)
log.addHandler(handler)
log.setLevel(logging.DEBUG)
log.propagate = False

log = logging.getLogger(__name__)

try:
    import tqdm
except:
    log.warning("For progress bar you need to pip install tqdm")
