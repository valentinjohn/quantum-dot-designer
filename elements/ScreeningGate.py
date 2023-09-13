# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from base import Element
import numpy as np
from helpers.helpers import create_segmented_path, get_polygons_from_path
import gdstk


class ScreeningGate(Element):
    def __init__(self, name):
        super().__init__(name)
        self.layer = 5
        self.vertices = []
        self.screen_paths = []
        self.qda_elements = None
        self._screen_path = False
        self.contact_vertices = None

    def screen(self, element_name, element_number,
               start_length, end_length,
               width=50e-3, relative_width=False):
        fo_line_name = f'fo_line_{element_name}_{element_number}'
        element_fo_path = self.qda_elements.components[fo_line_name].path
        screen_path = create_segmented_path(element_fo_path,
                                            start_length, end_length,
                                            width=width,
                                            relative_width=relative_width)
        self.screen_paths.append(screen_path)
        poly_vertices = get_polygons_from_path(screen_path)
        self.vertices.append(poly_vertices)

        self._screen_path = True

    def build(self):
        cell = gdstk.Cell(self.name)
        for vertices in self.vertices:
            screen = gdstk.Polygon(vertices,
                                   layer=self.layer)
            # screen.translate(self.x, self.y)
            screen.fillet(self.fillet, tolerance=self.fillet_tolerance)
            cell.add(screen)
            self.elements[self.name]['vertices'] = screen.points
            # self.elements[self.name]['positions'] = [[self.x, self.y]]
            self.elements[self.name]['positions'] = [[0, 0]]
            self.elements[self.name]['layer'] = self.layer

        self.cell = cell
        self._set_built(True)
