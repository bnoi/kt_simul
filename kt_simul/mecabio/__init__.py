"""
Mecabio module provides an easy API to model any bio-mechanical system.
"""

coords = ['x', 'y', 'z']
speed_coords = ['v'+c for c in coords]
dcoords = ['d'+c for c in coords]
ucoords = ['u'+c for c in coords]

# Import components
from .components import Structure
from .components import Point
from .components import Link

# Import model functions
from .dynamics import Model
from .dynamics import viscous
from .dynamics import dashpot
from .dynamics import spring
from .dynamics import dampedspring
from .dynamics import contraction
from .dynamics import linear_fv
