# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from QuantumDotDesigner.base import Element, Layer
from QuantumDotDesigner.BaseCollection import BaseCollection
import numpy as np
from QuantumDotDesigner.helpers.helpers import gen_poly
import gdstk
import copy


class BasicPolygon(Element):
    """
    Represents a basic polygon element in the QuantumDotDesigner system, extending the functionalities provided by the Element class.

    The BasicPolygon class is tailored for the creation of regular, standardized polygons. It allows for the definition of the number of corners, creating shapes ranging from triangles (minimal) to polygons with several sides. The class maintains the properties of standard elements, including layer association, scaling, and fillet operations.

    Unique to the BasicPolygon is its simple approach to creating regular shapes based on a single defining parameter: the number of corners. The build process generates a polygon with the specified number of corners, applying scaling and fillet operations to achieve the final geometric representation.

    Attributes:
        layer (Layer): The layer to which the basic polygon belongs. It must be a valid Layer instance.
        _layer_stage (str): The stage of the layer associated with the basic polygon, influencing its build process.
        diameter (float): The overall size attribute of the basic polygon, influencing its scale.
        corners (int): The number of corners (vertices) for the polygon, defining its basic shape.
        _asymx (float): Internal attribute to manage the asymmetry in the x-direction.
        _asymy (float): Internal attribute to manage the asymmetry in the y-direction, calculated based on _asymx.

    Methods:
        build(): Constructs the geometric representation of the basic polygon based on its attributes and adds it to a cell.
    """

    def __init__(self, name: str, collection: BaseCollection):

        super().__init__(name, collection)
        self.layer = None
        self._layer_stage = 'fine'
        self.diameter = None
        self.corners = 8
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
        Construct the basic polygon's geometric representation.

        This method generates a polygon based on the predefined number of corners, forming a regular shape. It applies scaling operations to adjust the polygon's size based on the diameter and asymmetry factors and performs fillet operations for edge smoothing. The final shape is then registered within a cell, and its geometric details are stored within the element's attributes.

        The build process is specific to the BasicPolygon class, enabling the straightforward creation of regular, standardized shapes.

        Raises:
            ValueError: If no valid Layer is assigned before building.
        """
        if self.layer is None or not isinstance(self.layer, Layer):
            raise ValueError(
                "A valid Layer must be assigned before building the element.")
        pl_points = gen_poly(self.corners)
        layer = getattr(self.layer, self._layer_stage)
        pl = gdstk.Polygon(pl_points, layer=layer)
        pl.scale(0.5 / np.cos(np.pi / 8) * self.diameter)
        pl.scale(sx=self._asymx, sy=self._asymy)
        pl.fillet(self.fillet, tolerance=self.fillet_tolerance)
        pl.fillet(0.02, tolerance=1e-4)
        cell = gdstk.Cell(self.name)
        cell.add(pl)
        self.elements[self.name]['vertices'] = pl.points
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.elements[self.name]['layer_stage'] = self._layer_stage
        self.cell = cell
        self._set_built(True)
