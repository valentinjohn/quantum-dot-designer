# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:41:19 2023

@author: vjohn
"""
# %% imports

from QuantumDotDesigner.base import PlotMixin
# from QuantumDotDesigner.elements.FanOutLineBase import FanOutLineBase
from QuantumDotDesigner.helpers.helpers import update_positions
import gdstk

# %% definition


class Sublattice(PlotMixin):
    def __init__(self, name):
        """
        Initialize a Sublattice object.

        Args:
            name (str): Name of the sublattice.
        """
        self.name = name
        self.elements = {}
        self.component = None
        self.rows = 1
        self.columns = 1
        self.spacing = (100, 100)
        self.center = (0, 0)
        self.cell = None
        self.xlim = None
        self.ylim = None
        self._width = (self.columns-1) * self.spacing[0]
        self._height = (self.rows-1) * self.spacing[1]
        self._built = False

    @property
    def built(self):
        return self._built

    # This setter is private and can only be used within the class
    def _set_built(self, new_value: bool):
        self._built = new_value

    def _update_width(self):
        """
        Update the width of the sublattice based on the number of columns and
        spacing.
        """
        self._width = (self.columns-1) * self.spacing[0]

    def _update_height(self):
        """
        Update the height of the sublattice based on the number of rows and
        spacing.
        """
        self._height = (self.rows-1) * self.spacing[1]

    def _get_lim(self, axis=0):
        """
        Get the minimum and maximum limits along the specified axis of the
        sublattice cell.

        Args:
            axis (int): Axis along which to find the limits (0 for x-axis,
                                                              1 for y-axis).

        Returns:
            tuple: Minimum and maximum limits as a tuple (min, max).
        """
        cell_max = 0
        cell_min = 0
        for poly in self.cell.get_polygons():
            poly_max = poly.points[:, axis].max()
            poly_min = poly.points[:, axis].min()
            cell_max = max(cell_max, poly_max)
            cell_min = min(cell_min, poly_min)
        if axis == 0:
            self.xlim = (cell_min, cell_max)
        elif axis == 1:
            self.ylim = (cell_min, cell_max)
        else:
            raise ValueError(
                f"axis has to be 0 or 1, but is '{axis}'.")

    def build(self):
        """
        Build the sublattice cell by adding the element references to the
        sublattice cell.
        """
        if not self.component.built:
            self.component.build()

        self._update_height()
        self._update_width()
        cell = gdstk.Cell(self.name)
        cell.add(gdstk.Reference(self.component.cell,
                                 (self.center[0] - self._width/2,
                                  self.center[1] - self._height/2),
                                 columns=self.columns,
                                 rows=self.rows,
                                 spacing=(self.spacing[0], self.spacing[1])
                                 ))
        self.cell = cell
        self._set_built(True)
        self._get_lim(axis=0)
        self._get_lim(axis=1)
        # only update elements attribute if it exists
        if hasattr(self, 'elements'):
            self.elements = update_positions(self.component.elements,
                                             self.rows, self.columns,
                                             self.spacing[0], self.spacing[1],
                                             self.center)
