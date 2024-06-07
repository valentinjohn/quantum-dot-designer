"""
Quantum Dot Designer is a Python package for designing quantum dot arrays. It is built on top of the GDSII library and provides a high-level API for designing quantum dot arrays. The package is designed to be modular and extensible, allowing users to easily create custom elements and components.
"""

__version__ = "1.1.2"

from .QuantumDotArray import QuantumDotArray
from .BaseCollection import BaseCollection

from . import elements
from . import components


from .fanout import Fanout, FanoutPoints
from .base.UnitCell import UnitCell

