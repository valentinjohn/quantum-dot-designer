# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:37:28 2023

@author: vjohn
"""
# %% imports

from base import PlotMixin
from base.Sublattice import Sublattice
from helpers.helpers import merge_device_positions
import gdstk

# %% definition


class UnitCell(PlotMixin):
    def __init__(self, name):
        """
        Initialize a UnitCell object.

        Args:
            parent_instance: The QuantumDotArray instance that this UnitCell belongs to.
            name (str): Name of the unit cell (default: 'unit_cell').
        """

        self.name = name
        self.elements = {}
        self.components = {}
        self.cell = gdstk.Cell(name)
        self._xlim = None
        self._ylim = None
        self._n = 0
        self._built = False

    @property
    def built(self):
        return self._built

    # This setter is private and can only be used within the class
    def _set_built(self, new_value: bool):
        self._built = new_value

    def add_component(self, component=None, build=False):
        """
        Add a sublattice to the unit cell.

        Args:
            name (str): Name of the sublattice.
            component: You can assign otionally a component in the argument
            build: Adds and builds at the same time

        Returns:
            Sublattice: The created Sublattice object.
        """
        name = f'{self.name}_sublattice_{self._n}'
        self._n = self._n + 1
        sublattice = Sublattice(name)
        if component is not None:
            if component._built == False:
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

    def _get_lim(self, axis=0):
        cell_max = 0
        cell_min = 0
        for poly in self.cell.get_polygons():
            poly_max = poly.points[:, axis].max()
            poly_min = poly.points[:, axis].min()
            cell_max = max(cell_max, poly_max)
            cell_min = min(cell_min, poly_min)
        if axis:
            self._ylim = (cell_min, cell_max)
        else:
            self._xlim = (cell_min, cell_max)

    @property
    def xlim(self):
        if not self.built:
            self.build()
        return self._xlim

    @property
    def ylim(self):
        if not self.built:
            self.build()
        return self._ylim

    def build(self):
        """
        Build the unit cell by adding sublattices to the unit cell.

        This method adds the sublattices to the unit cell's cell object.
        """

        elements = {}
        for cell in self.components.values():
            if not cell._built:
                cell.build()
            self.cell.add(gdstk.Reference(cell.cell))
            elements = merge_device_positions(elements, cell.elements)
        self.elements = elements
        self.cell.flatten()
        self._get_lim(axis=0)
        self._get_lim(axis=1)
        self._set_built(True)
