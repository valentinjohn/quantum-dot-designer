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


class Plunger(Element):
    """
    Represents a plunger element in the QuantumDotDesigner system, extending the functionalities provided by the Element class.

    The Plunger class is designed for elements that require specific geometric constructions, particularly a circular shape with adjustable asymmetry in the x and y directions. It is associated with a specific layer and has a defined diameter that influences its build process.

    Unique to the Plunger is its ability to adjust its asymmetry properties, allowing for non-uniform scaling in its geometric representation. The build process for the Plunger involves creating a polygon with specific attributes and adjustments, followed by scaling and fillet operations to achieve the desired final shape.

    Attributes:
        layer (Layer): Identifier for the layer to which the plunger belongs. Can be None if the plunger is not associated with a layer.
        layer_stage (str): The stage of the layer associated with the plunger, influencing its build process.
        diameter (type): The overall size attribute of the plunger, influencing its scale.
        _asymx (float): Private attribute to manage the asymmetry in the x-direction.
        _asymy (float): Private attribute to manage the asymmetry in the y-direction, calculated based on _asymx.
        corners (int): number of corners of the polygon. Defaults to 8.

    Methods:
        build(): Constructs the geometric representation of the plunger based on its attributes and adds it to a cell.
    """

    def __init__(self, name: str, collection: BaseCollection):
        super().__init__(name, collection)
        self.layer = None
        self._layer_stage = 'fine'
        self.diameter = None
        self._asymx = 1
        self._asymy = 1 / self._asymx
        self.corners = 8

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
        Construct the plunger's geometric representation.

        This method generates a polygon that forms the basis of the plunger, applies scaling and fillet operations 
        to adjust its shape and size, and registers the final geometry and its properties to the associated cell. 
        The method relies on the plunger's attributes such as diameter, layer information, and asymmetry factors 
        to create the intended design.

        The build process is specific to the Plunger class and takes into account its unique requirements and properties.

        Raises:
            ValueError: If no valid Layer is associated with the plunger before building.
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
