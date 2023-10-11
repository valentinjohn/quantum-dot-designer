# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from QuantumDotDesigner.base import Element
from QuantumDotDesigner.BaseCollection import BaseCollection
import numpy as np
import gdstk


class ArbitraryPolygon(Element):
    def __init__(self, name: str, collection: BaseCollection):
        """
        Initialize an ArbitraryPolygon object.

        Args:
            name (str): Name of the element.
            layer (int): Description of layer. Default is None.
            rotate (float): Rotation of the polyunger. Default is 0.0.
            fillet (float): Fillet value. Indicates how much the polygon is rounded.
            ... (other attributes)
        """
        super().__init__(name, collection)
        self.layer = None
        self.layer_stage = 'fine'
        self.vertices = None
        self._asymx = 1
        self._asymy = 1 / self._asymx

    @property
    def asymx(self):
        return self._asymx

    @asymx.setter
    def asymx(self, value):
        if value != 0:  # to avoid division by zero
            self._asymx = value
            self._asymy = 1 / self._asymx

    @property
    def asymy(self):
        return self._asymy

    @asymy.setter
    def asymy(self, value):
        if value != 0:  # to avoid division by zero
            self._asymy = value
            self._asymx = 1 / self._asymy

    def build(self):
        """
        Build the polyunger element.
        """
        layer = getattr(self.layer, self.layer_stage)
        poly = gdstk.Polygon(self.vertices, layer=layer)
        poly.scale(0.5 / np.cos(np.pi / 8) * self.diameter)
        poly.scale(sx=self._asymx, sy=self._asymy)
        poly.fillet(self.fillet, tolerance=self.fillet_tolerance)
        poly.fillet(0.02, tolerance=1e-4)
        cell = gdstk.Cell(self.name)
        cell.add(poly)
        self.elements[self.name]['vertices'] = poly.points
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.elements[self.name]['layer_stage'] = self.layer_stage
        self.cell = cell
        self._set_built(True)
