# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""
from ..base import Element, Layer
from ..BaseCollection import BaseCollection


import gdstk
import copy


class Barrier(Element):
    """
    Represents a barrier element in the quantum_dot_designer system, extending the functionalities provided by the Element class.

    The Barrier class is tailored for elements that simulate physical barriers, having specific geometric attributes such as width and length. These barriers are rectangular with additional geometric features, and they can be associated with a specific layer.

    The construction of a barrier involves creating a polygon that represents the barrier's physical shape, including its orientation (rotation) and fillet, if any. This polygon is then registered within a cell, and its geometric details are stored within the element's attributes.

    Attributes:
        layer (Layer): Identifier for the layer to which the barrier belongs. Can be None if the barrier is not associated with a layer.
        _layer_stage (str): Private attribute indicating the stage of the layer associated with the barrier. Defaults to 'fine'.
        width (type): The horizontal dimension of the barrier.
        length (type): The vertical dimension of the barrier.

    Methods:
        build(): Constructs the geometric representation of the barrier based on its attributes and adds it to a cell.
    """

    def __init__(self, name, collection: BaseCollection):
        super().__init__(name, collection)
        self.layer = None
        self._layer_stage = 'fine'
        self.width = None
        self.length = None

    def build(self):
        """
        Construct the barrier's geometric representation.

        This method generates a polygon that forms the basis of the barrier, taking into account its width and length to create the desired shape. The barrier is constructed with certain distinctive vertices to represent its physical characteristics. The method then applies any specified rotation and fillet, finalizing the barrier's geometry.

        The constructed barrier is then added to a new cell, and its geometric and layer details are registered within the element's attributes. The build process is specific to the Barrier class and utilizes its unique properties.

        Raises:
            ValueError: If no valid Layer is associated with the barrier before building.
        """
        if self.layer is None or not isinstance(self.layer, Layer):
            raise ValueError(
                "A valid Layer must be assigned before building the element.")

        layer = getattr(self.layer, self._layer_stage)
        bar = gdstk.Polygon([(-self.length/2, -self.width/2),
                             (self.length/2, -self.width/2),
                             (self.length/2 + self.width/4, -self.width/8),
                             (self.length/2 + self.width/4, self.width/8),
                             (self.length/2, self.width/2),
                             (-self.length/2, self.width/2),
                             (-self.length/2 - self.width/2, self.width/8),
                             (-self.length/2 - self.width/2, -self.width/8)],
                            layer=layer)
        bar.rotate(self.rotate)
        bar.fillet(self.fillet, tolerance=self.fillet_tolerance)
        cell = gdstk.Cell(self.name)
        cell.add(bar)
        self.elements[self.name]['vertices'] = bar.points
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.elements[self.name]['layer_stage'] = self._layer_stage
        self.cell = cell
        self._set_built(True)
