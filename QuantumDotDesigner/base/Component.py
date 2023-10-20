# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 09:37:32 2023

@author: vjohn
"""

# %% imports

from QuantumDotDesigner.base import UnitCell

# %% definition


class Component(UnitCell):
    """
    Represents a component within a quantum dot design, extending the functionality of a unit cell.

    The `Component` class is derived from the `UnitCell` class and is intended to encapsulate a set of elements and/or other components, providing additional layers of organization or functionality within a quantum dot array. While currently serving as a direct pass-through to `UnitCell`, it exists as a separate class for future expansion and clearer semantic representation in a quantum dot system.

    Attributes inherited from `UnitCell`:
        name (str): Unique identifier for the component/unit cell.
        elements (dict): Dictionary of elements that constitute the component/unit cell.
        components (dict): Dictionary of sub-components or sublattices added to the component/unit cell.
        cell (gdstk.Cell): Underlying representation of the component/unit cell in GDSTK.
        _xlim (tuple, optional): Horizontal (x-axis) limits of the component/unit cell's geometry.
        _ylim (tuple, optional): Vertical (y-axis) limits of the component/unit cell's geometry.
        _n (int): Counter used for naming sublattices.
        _built (bool): Flag indicating if the component/unit cell has been constructed.
    """
    pass
