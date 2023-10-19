# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from QuantumDotDesigner.base import Element, Layer
from QuantumDotDesigner.BaseCollection import BaseCollection
import numpy as np
import gdstk


class ArbitraryPolygon(Element):
    """
    Represents an arbitrary polygon element in the QuantumDotDesigner system, extending the functionalities provided by the Element class.

    The ArbitraryPolygon class caters to elements that require custom geometric shapes, allowing users to define a polygon with an arbitrary number of vertices. This class maintains the properties of standard elements while providing additional flexibility in terms of shape design.

    Unique to the ArbitraryPolygon is its reliance on a set of vertices for its shape definition, rather than predetermined geometric parameters. It also allows for non-uniform scaling via asymmetry properties. The build process involves creating a polygon based on the provided vertices, followed by scaling, rotation, and fillet operations to achieve the desired final shape.

    Attributes:
        layer (Layer): The layer to which the arbitrary polygon belongs. It must be a valid Layer instance.
        _layer_stage (str): The stage of the layer associated with the arbitrary polygon, defaults to 'fine'.
        vertices (list): A list of (x, y) tuples that define the vertices of the polygon. These vertices form the perimeter of the arbitrary shape.
        _asymx (float): Internal attribute to manage the asymmetry in the x-direction.
        _asymy (float): Internal attribute to manage the asymmetry in the y-direction, calculated based on _asymx.

    Methods:
        build(): Constructs the geometric representation of the arbitrary polygon based on its vertices and adds it to a cell.
    """

    def __init__(self, name: str, collection: BaseCollection):
        super().__init__(name, collection)
        self.layer = None
        self._layer_stage = 'fine'
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
        Construct the arbitrary polygon's geometric representation.

        This method takes the predefined set of vertices and forms a polygon. It applies scaling operations to adjust the polygon's shape based on the asymmetry factors and performs fillet operations for edge smoothing. The final shape is then registered within a cell, and its geometric details are stored within the element's attributes.

        The build process is specific to the ArbitraryPolygon class and allows for the creation of custom shapes beyond standard geometric figures.

        Raises:
            ValueError: If no valid Layer is assigned before building.
        """
        if self.layer is None or not isinstance(self.layer, Layer):
            raise ValueError(
                "A valid Layer must be assigned before building the element.")
        layer = getattr(self.layer, self._layer_stage)
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
        self.elements[self.name]['layer_stage'] = self._layer_stage
        self.cell = cell
        self._set_built(True)
