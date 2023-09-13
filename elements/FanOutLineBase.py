# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from base import PlotMixin
import gdstk


class FanOutLineBase(PlotMixin):
    def __init__(self, name):
        self.name = name
        self.element_name = None
        self.element_number = None
        self.fo_fine_coarse_overlap = None
        self.fo_fine_coarse_overlap_gap = 0.3
        self.layer = None
        self.polygons = None
        self.path = None
        self.fo_direction = None
        self.n_fanout = None
        self.cell = gdstk.Cell(self.name)
        self.fillet = 0
        self.fillet_tolerance = 1e-3
        self.elements = {self.name: {'vertices': [],
                                     'positions': [],
                                     'layer': self.layer}}
        self._built = False

    @property
    def built(self):
        return self._built

    # This setter is private and can only be used within the class
    def _set_built(self, new_value: bool):
        self._built = new_value

    def build(self):
        fo_line = gdstk.Polygon(self.polygons, layer=self.layer)
        fo_line.fillet(self.fillet, tolerance=self.fillet_tolerance)
        fo_line.fillet(0.02, tolerance=1e-4)

        self.elements[self.name]['vertices'] = fo_line.points
        self.elements[self.name]['positions'] = [0, 0]
        self.elements[self.name]['layer'] = self.layer

        self.cell.add(fo_line)
        self._set_built(True)
