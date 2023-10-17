# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:07:21 2023

@author: vjohn
"""

# %% imports

from abc import abstractmethod
from QuantumDotDesigner.base import PlotMixin

# %% definition


class ElementBase(PlotMixin):
    def __init__(self, name: str):
        """
        Initialize an Element object.

        Args:
            name (str): Name of the element.
            layer (int): Description of layer. Default is None.
            rotate (float): Rotation of the element. Default is 0.0.
            fillet (float): Fillet value. Indicates how much the polygon is rounded.
        """

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
