from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

def in_ipython():
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True
