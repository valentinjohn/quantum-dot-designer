# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:07:21 2023

@author: vjohn
"""

# %% imports

from abc import abstractmethod
from .PlotMixin import PlotMixin

# %% definition


class ElementBase(PlotMixin):
    """
    Abstract base class for elements in the quantum_dot_designer system, providing a common interface,
    properties, and methods for all element objects. This class should not be instantiated directly 
    but should be subclassed to create specific element types with their own implementations.

    Each element has a unique name and contains information about its geometry, layer, and additional 
    properties. Elements can be manipulated through various parameters including rotation and fillet, 
    and they need to be constructed using the `build` method implemented specifically for each subtype.

    The class also supports creating a copy of an element with the possibility to override specific attributes.

    Attributes:
        name (str): Unique identifier for the element.
        layer (Layer): Layer object associated to the element. Defaults to None.
        elements (dict): Dictionary containing geometric information about the built elements.
        cell (gdstk.cell, optional): gdstk cell. Defaults to None.
        _fo_contact_point (tuple): The fanout contact point in the form (x, y).
        _fo_contact_width (float): The width of the fanout contact point.
        _fo_contact_vector (vector, optional): Unit vector indicating the direction, into which the contact has to fanout.
        rotate (float): Rotation angle of the element in degrees. Defaults to 0.0.
        fillet (float): Degree to which the element's polygon corners are rounded. Defaults to 0.0.
        fillet_tolerance (float): Tolerance applied during the fillet operation.
        _built (bool): Flag indicating whether the element has been built. Defaults to False.
        bondpad_off (bool): Flag indicating whether the bondpads have to be placed off the regular bondpad strip. Defaults to False.

    Methods:
        copy(copy_name: str): Create a copy of the element with a new name, maintaining the original element's properties.
        build(): Abstract method. Constructs the element's geometry. Must be implemented by subclasses.
    """

    def __init__(self, name: str):
        self.name = name
        self.layer = None
        self.elements = {name: {'vertices': [],
                                'positions': [],
                                'layer': self.layer}}
        self.cell = None
        self._fo_contact_point = (0, 0)
        self._fo_contact_width = 40e-3
        self._fo_contact_vector = None
        self.rotate = 0.0
        self.fillet = 0.0
        self.fillet_tolerance = 1e-3
        self._built = False
        self.bondpad_off = False

    # default attributes to skip for Element
    _skip_copy_attrs = {'name', 'cell', 'elements'}

    @property
    def built(self):
        return self._built

    # This setter is private and can only be used within the class
    def _set_built(self, new_value: bool):
        self._built = new_value

    def copy(self, copy_name):
        # Use the same class as the current instance
        copied_obj = self.__class__(copy_name)

        for attr, value in self.__dict__.items():
            if attr not in self._skip_copy_attrs:
                setattr(copied_obj, attr, value)

        return copied_obj

    @abstractmethod
    def build(self):
        """
        Build the element.

        This method should be implemented in subclasses to build the specific element.
        """
        pass
