# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from ..base import Element, Layer
from ..BaseCollection import BaseCollection
import numpy as np
from ..helpers.helpers import (rot_mat, midpoint,
                                                  orthogonal_unit_vector,
                                                  distance,
                                                  adjust_vector_direction)

import gdstk


class Ohmic(Element):
    """
    Represents an ohmic contact within the quantum_dot_designer system.

    It is an extension of the standard Element class, equipped with additional attributes specific to ohmic contacts, such as contact length, contact offset, and contact angle.

    Attributes:
        layer (Layer): The layer to which the ohmic contact belongs. It must be a valid Layer instance.
        _layer_stage (str): The stage of the layer associated with the ohmic contact, defaults to 'fine'.
        contact_length (float): The length of the ohmic contact.
        contact_offset (float): The offset distance of the ohmic contact from a reference point.
        contact_angle (float): The angle at which the ohmic contact is oriented.
        contact_width (float): The width of the ohmic contact.
        sensor_pos (str): Position of the sensor relative to the device.
        ohmic_pos (str): Position of the ohmic contact relative to the sensor.
        rotate (float): The angle of rotation applied to the ohmic contact.
        vertices (list): A list of (x, y) tuples defining the vertices of the ohmic contact.
        bondpad_off (bool): A flag indicating whether the bondpad is off. Defaults to True.

    Methods:
        compute_ohmic_vertices(): Calculates the vertices of the ohmic contact based on geometric properties.
        build(): Constructs the geometric representation of the ohmic contact and adds it to a cell.
    """

    def __init__(self, name, collection: BaseCollection):
        super().__init__(name, collection)
        self.layer = None
        self._layer_stage = 'fine'
        self.contact_length = 0.095
        self.contact_offset = 0.01
        self.contact_angle = np.pi/4
        self.contact_width = 0.05
        self.sensor_pos = 'top'
        self.ohmic_pos = 'right'
        self.rotate = 0
        self.vertices = None
        self.bondpad_off = True

    def compute_ohmic_vertices(self):
        """
        Calculate the vertices defining the ohmic contact.

        This method uses the geometric properties of the ohmic contact, such as its length, width, angle, and position, 
        to calculate the vertices that define its shape. These vertices are essential for constructing the ohmic contact's 
        accurate geometric representation.
        """
        multiplier_dict = {('top', 'right'): 1,
                           ('top', 'left'): -1,
                           ('bottom', 'right'): -1,
                           ('bottom', 'left'): 1,
                           ('right', 'top'): -1,
                           ('right', 'bottom'): 1,
                           ('left', 'top'): 1,
                           ('left', 'bottom'): -1
                           }

        rotation_dict = {'top': 0,
                         'bottom': np.pi,
                         'right': -np.pi/2,
                         'left': np.pi/2,
                         }

        multiplier = multiplier_dict[(self.sensor_pos, self.ohmic_pos)]
        rotation = rotation_dict[self.sensor_pos]

        v1 = [0, self.contact_length/2]
        v2 = (0, -self.contact_length/2)
        v4 = (multiplier * self.contact_offset, self.contact_length/2)
        v3 = np.array(v4) + self.contact_width * \
            np.array([multiplier * np.cos(self.contact_angle), -
                     np.sin(self.contact_angle)])

        v1 = list(rot_mat(rotation) @ np.array(v1))
        v2 = list(rot_mat(rotation) @ np.array(v2))
        v3 = list(rot_mat(rotation) @ np.array(v3))
        v4 = list(rot_mat(rotation) @ np.array(v4))

        vertices = [v1, v2, v3, v4]

        self.vertices = vertices

    def build(self):
        """
        Construct the ohmic contact's geometric representation.

        This method generates a polygon representing the ohmic contact based on the calculated vertices. It applies the 
        necessary geometric transformations, including rotation and fillet, to finalize the shape. The method then registers 
        the ohmic contact within a cell.

        The build process is specific to the Ohmic class and takes into account the contact's unique geometric requirements 
        and properties.

        Raises:
            ValueError: If no valid Layer is assigned before building.
        """
        if self.layer is None or not isinstance(self.layer, Layer):
            raise ValueError(
                "A valid Layer must be assigned before building the element.")

        self.compute_ohmic_vertices()
        layer = getattr(self.layer, self._layer_stage)
        ohmic = gdstk.Polygon(self.vertices, layer=layer)
        ohmic.fillet(self.fillet, tolerance=self.fillet_tolerance)
        ohmic.rotate(self.rotate)
        self._fo_contact_point = midpoint(ohmic.points[-2],
                                          ohmic.points[-1])
        vector = orthogonal_unit_vector(ohmic.points[-2],
                                        ohmic.points[-1])
        point = -ohmic.points[-3]
        self._fo_contact_vector = adjust_vector_direction(vector, point)
        self._fo_contact_width = distance(ohmic.points[-2],
                                          ohmic.points[-1])
        cell = gdstk.Cell(self.name)
        cell.add(ohmic)
        self.elements[self.name]['vertices'] = ohmic.points
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.elements[self.name]['layer_stage'] = self._layer_stage
        self.cell = cell
        self._set_built(True)
