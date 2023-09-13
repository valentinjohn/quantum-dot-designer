# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from QuantumDotDesigner.base import Element
import gdstk


class Barrier(Element):
    def __init__(self, name):
        """
        Initialize a Barrier object.

        Args:
            name (str): Name of the barrier.
        """
        super().__init__(name)
        self.layer = 5
        self.width = None
        self.length = None

    def build(self):
        """
        Build the barrier element.
        """
        bar = gdstk.Polygon([(-self.length/2, -self.width/2),
                             (self.length/2, -self.width/2),
                             (self.length/2 + self.width/4, -self.width/8),
                             (self.length/2 + self.width/4, self.width/8),
                             (self.length/2, self.width/2),
                             (-self.length/2, self.width/2),
                             (-self.length/2 - self.width/2, self.width/8),
                             (-self.length/2 - self.width/2, -self.width/8)],
                            layer=self.layer)
        bar.rotate(self.rotate)
        # bar.translate(self.x, self.y)
        bar.fillet(self.fillet, tolerance=self.fillet_tolerance)
        cell = gdstk.Cell(self.name)
        cell.add(bar)
        self.elements[self.name]['vertices'] = bar.points
        # self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.cell = cell
        self._set_built(True)
