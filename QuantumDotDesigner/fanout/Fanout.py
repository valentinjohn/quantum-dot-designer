# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 10:11:40 2023

@author: vjohn
"""

import gdstk

from QuantumDotDesigner.helpers.helpers import merge_device_positions
from QuantumDotDesigner.base import Sublattice, PlotMixin


class Fanout(PlotMixin):
    def __init__(self, name):
        self.name = name
        self.elements = {}
        self.cell = gdstk.Cell(name)
        self.components = {}
        self._built = False
        self._n = 0

    @property
    def built(self):
        return self._built

    def add_component(self, component=None, build=False):
        """
        Add a sublattice to the unit cell.

        Args:
            name (str): Name of the sublattice.
            component: You can assign optionally a component in the argument
            build: Adds and builds at the same time with a single component placed at origin.

        Returns:
            Sublattice: The created Sublattice object.
        """
        name = f'{self.name}_sublattice_{self._n}'
        self._n = self._n + 1
        sublattice = Sublattice(name)
        if component is not None:
            if not component._built:
                component.build()
            sublattice.component = component
            if build:
                sublattice.build()
            else:
                pass
        else:
            pass
        self.components[name] = sublattice
        return sublattice

    def build(self):
        elements = {}
        for cell in self.components.values():
            if not cell._built:
                cell.build()
            self.cell.add(gdstk.Reference(cell.component.cell))
            elements = merge_device_positions(
                elements, cell.component.elements)
        self.elements = elements
        self.cell.flatten()
        self._built = True
