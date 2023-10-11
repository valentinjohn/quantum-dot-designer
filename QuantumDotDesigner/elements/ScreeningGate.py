# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from QuantumDotDesigner.base import Element
from QuantumDotDesigner.BaseCollection import BaseCollection
import numpy as np
from QuantumDotDesigner.helpers.helpers import (create_segmented_path,
                                                get_polygons_from_path,
                                                get_end_path,
                                                point_along_line,
                                                orthogonal_unit_vector,
                                                adjust_vector_direction)
import gdstk


class ScreeningGate(Element):
    def __init__(self, name, collection: BaseCollection):
        super().__init__(name, collection)
        self.layer = None
        self.layer_stage = 'fine'
        self.vertices = []
        self.screen_paths = []
        self.collection = collection
        self._screen_path = False
        self.contact_vertices = None
        self.fo_contact_width = 40e-3
        self.fo_contact_direction = 0

    def screen(self, element_name, element_number,
               points, widths):
        fo_line_name = f'fo_line_{element_name}_{element_number}'
        element_fo_path = self.collection.elements[fo_line_name].path
        screen_path = create_segmented_path(element_fo_path,
                                            points,
                                            widths)
        self.screen_paths.append(screen_path)
        poly_vertices = get_polygons_from_path(screen_path)
        self.vertices.append(poly_vertices)

        self._screen_path = True

    def build(self):
        cell = gdstk.Cell(self.name)
        layer = getattr(self.layer, self.layer_stage)
        for vertices in self.vertices:
            screen = gdstk.Polygon(vertices,
                                   layer=layer)
            # screen.translate(self.x, self.y)
            screen.fillet(self.fillet, tolerance=self.fillet_tolerance)
            cell.add(screen)
            self.elements[self.name]['vertices'] = screen.points
            # self.elements[self.name]['positions'] = [[self.x, self.y]]
            self.elements[self.name]['positions'] = [[0, 0]]
            self.elements[self.name]['layer'] = self.layer

        if self.screen_paths:
            end_path = get_end_path(self.vertices)
            self._fo_contact_point = point_along_line(end_path[self.fo_contact_direction % 2],
                                                      end_path[(
                                                          self.fo_contact_direction+1) % 2],
                                                      self.fo_contact_width/2)
            vector = orthogonal_unit_vector(end_path[0],
                                            end_path[1])
            point = self.screen_paths[0][0][:2]
            self._fo_contact_vector = adjust_vector_direction(vector, point)

        self.cell = cell
        self._set_built(True)
