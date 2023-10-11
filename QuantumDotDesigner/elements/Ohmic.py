# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from QuantumDotDesigner.base import Element
from QuantumDotDesigner.BaseCollection import BaseCollection
import numpy as np
from QuantumDotDesigner.helpers.helpers import (rot_mat, midpoint,
                                                orthogonal_unit_vector,
                                                distance,
                                                adjust_vector_direction)
import gdstk


class Ohmic(Element):
    def __init__(self, name, collection: BaseCollection):
        super().__init__(name, collection)
        self.layer = None
        self.layer_stage = 'fine'
        self.contact_length = 0.095
        self.contact_offset = 0.01
        self.contact_angle = np.pi/4
        self.contact_width = 0.05
        self.sensor_pos = 'top'
        self.ohmic_pos = 'right'
        self.rotate = 0
        self.vertices = None

    def compute_ohmic_vertices(self):
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
        Build the ohmic element.
        """
        self.compute_ohmic_vertices()
        layer = getattr(self.layer, self.layer_stage)
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
        # self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.elements[self.name]['layer_stage'] = self.layer_stage
        self.cell = cell
        self._set_built(True)
