# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:53:24 2023

@author: vjohn
"""
# %% imports

from QuantumDotDesigner.base import UnitCell
from QuantumDotDesigner.QuantumDotArrayElements import QuantumDotArrayElements
import gdstk

# %% definition


class FanOutLine(UnitCell):
    def __init__(self, element_name: str, element_number: int,
                 qda_elements: QuantumDotArrayElements):
        super().__init__(f'fo_{element_name}_{element_number}')
        self.qda_elements = qda_elements
        self.fo_points = None
        name_coarse = f'fo_coarse_{element_name}_{element_number}'
        self.fo_line_coarse = qda_elements.add_fo_line_coarse(name_coarse)
        name_fine = f'fo_line_{element_name}_{element_number}'
        self.fo_line_fine = qda_elements.add_fo_line_fine(name_fine)
        self.element_name = element_name
        self.element_number = element_number
        # self.fo_fine_coarse_overlap = None
        # self.fo_fine_coarse_overlap_gap = 0.3
        self.element = qda_elements.components[element_name]
        self.start_offset = self.element.fo_contact_point
        self.fo_line_fine.fo_width_start = self.element.fo_contact_width
        self.fo_direction_start = self.element.fo_contact_vector
        self.layer = self.element.layer
        self.polygons = None
        self.fo_direction = None
        self.n_fanout = None
        self.cell = gdstk.Cell(self.name)
        self.fillet = 0
        self.fillet_tolerance = 1e-3
        self.elements = {}
        # self.fine_fo_width_start = 40e-3
        # self.fine_fo_start = None
        # self.fine_fo_end = None
        # self.points_along_path = []
        self._built = False

        self.fo_line_fine.element_name = self.element_name
        self.fo_line_fine.element_number = self.element_number

        self.fo_line_coarse.element_name = self.element_name
        self.fo_line_coarse.element_number = self.element_number

    # def update_properties(self):
    #     self.element_name =

    def add_coarse_fo_line(self):

        # self.fo_line_coarse.name = name
        # self.fo_line_coarse.cell.name = name
        self.elements[self.fo_line_coarse.name] = {'vertices': [],
                                                   'positions': [],
                                                   'layer': self.fo_line_coarse.layer}
        # self.fo_line_coarse = self.qda_elements.add_fo_line_coarse(name)

        self.fo_line_coarse.fo_direction = self.fo_direction
        self.fo_line_coarse.layer = self.qda_elements.components[self.element_name].layer + 20

        attributes = {'fo_direction': self.fo_direction,
                      'n_fanout': self.n_fanout,
                      'fo_points': self.fo_points}
        missing_attrs = [name for name,
                         value in attributes.items() if value is None]
        if missing_attrs:
            raise ValueError(
                f"Attributes {', '.join(missing_attrs)} cannot be None.")
        else:
            self.fo_line_coarse.polygons = self.fo_points.fo_polygons_coarse[
                self.fo_direction][self.n_fanout]

        # self.components[name] = self.fo_line_coarse
        self.add_component(self.fo_line_coarse, build=True)

    def add_fine_fo_line(self):
        # name = f'fo_line_{self.element_name}_{self.element_number}_fine'
        # self.fo_line_fine.name = name
        # self.fo_line_fine.cell.name = name
        if self.fo_direction_start:
            self.fo_line_fine.points_along_path = [[0.01*self.fo_direction_start[0],
                                                   0.01 *
                                                   self.fo_direction_start[1],
                                                   'start']]

        self.elements[self.fo_line_fine.name] = {'vertices': [],
                                                 'positions': [],
                                                 'layer': self.fo_line_fine.layer}

        self.fo_line_fine.fo_direction = self.fo_direction
        self.fo_line_fine.layer = self.qda_elements.components[self.element_name].layer
        self.fo_line_fine.fillet = self.fillet

        attributes = {
            'fo_direction': self.fo_direction,
            'n_fanout': self.n_fanout,
            'fo_points': self.fo_points
        }

        missing_attrs = [name for name,
                         value in attributes.items() if value is None]

        if missing_attrs:
            raise ValueError(
                f"Attributes {', '.join(missing_attrs)} cannot be None.")
        else:
            self.fo_line_fine.fo_start = self.fo_points.qda.elements[
                self.element_name]['positions'][self.element_number]
            self.fo_line_fine.fo_start[0] += self.start_offset[0]
            self.fo_line_fine.fo_start[1] += self.start_offset[1]
            self.fo_line_fine.fo_end = self.fo_points.get_fo_overlap_points(
                self.n_fanout, self.fo_direction)

        # self.components[name] = self.fo_line_fine
        self.add_component(self.fo_line_fine, build=True)

    def build(self):
        self.add_coarse_fo_line()
        self.add_fine_fo_line()
        super().build()
