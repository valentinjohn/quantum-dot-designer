# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from QuantumDotDesigner.base import Element
from QuantumDotDesigner.helpers.helpers import generate_clavier_gates
from QuantumDotDesigner.BaseCollection import BaseCollection
import gdstk


class ClavierGate(Element):
    def __init__(self, name, collection: BaseCollection):
        """
        Initialize a Clavier gate object.

        Args:
            name (str): Name of the clavier gate
        """
        super().__init__(name, collection)
        self.layer = None
        self.layer_stage = 'fine'
        self.width = 100
        self.length = 100*4*10
        self.gate_width = 100
        self.gate_length = 300
        self.n_clav_rep = 10
        self.spacing = 100*4
        self.shift = 0
        self.x = 0
        self.y = 0
        self.fillet = 0  # 0.02
        self.fillet_tolerance = 1e-4
        self.rotation = 0

    def build(self):
        """
        Build the clavier gate element.
        """
        cl_points = generate_clavier_gates(self.width, self.length,
                                           self.gate_width,
                                           self.gate_length,
                                           self.n_clav_rep, self.spacing,
                                           self.shift, (-self.length/2, 0),
                                           self.rotation)
        layer = getattr(self.layer, self.layer_stage)
        cl = gdstk.Polygon(cl_points, layer=layer)

        cl.translate(self.x, self.y)
        cl.fillet(self.fillet, tolerance=self.fillet_tolerance)
        cell = gdstk.Cell(self.name)
        cell.add(cl)
        self.elements[self.name]['vertices'] = cl.points
        self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['layer'] = self.layer
        self.cell = cell
        self._set_built(True)
