# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:53:24 2023

@author: vjohn
"""
# %% imports

from QuantumDotDesigner.base import Component
from QuantumDotDesigner.BaseCollection import BaseCollection
from QuantumDotDesigner.fanout import FanoutPoints

from QuantumDotDesigner.elements.Ohmic import Ohmic
from QuantumDotDesigner.elements.FanOutLineFine import FanOutLineFine
from QuantumDotDesigner.elements.FanOutLineCoarse import FanOutLineCoarse

import gdstk
import matplotlib.pylab as plt

# %% definition


class FanOutLine(Component):
    def __init__(self, element_name: str, element_number: int,
                 collection: BaseCollection, fo_points: FanoutPoints):
        # if not isinstance(qda_elements, QuantumDotArrayElements):
        #     raise TypeError(
        #         f"Expected qda_elements to be of type {QuantumDotArrayElements}, but got {type(qda_elements)} instead.")
        super().__init__(f'fo_{element_name}_{element_number}')
        self._init_elements(element_name, element_number, collection)
        self.collection = collection
        self.fo_points = fo_points
        self.element_name = element_name
        self.element_number = element_number
        # self.fo_fine_coarse_overlap = None
        # self.fo_fine_coarse_overlap_gap = 0.3
        self.element = self.collection.elements[element_name]
        self.start_offset = self.element._fo_contact_point
        self.fo_line_fine.fo_width_start = self.element._fo_contact_width
        self.fo_direction_start = self.element._fo_contact_vector
        self.layer = self.element.layer
        self.polygons = None
        self.fo_direction = None
        self.n_fanout = None
        self.cell = gdstk.Cell(self.name)
        self.fillet = 0
        self.fillet_tolerance = 1e-3
        # self.elements = {}
        # self.fine_fo_width_start = 40e-3
        # self.fine_fo_start = None
        # self.fine_fo_end = None
        # self.points_along_path = []
        self._built = False

        self.fo_line_fine.element_name = self.element_name
        self.fo_line_fine.element_number = self.element_number

        self.fo_line_coarse.element_name = self.element_name
        self.fo_line_coarse.element_number = self.element_number

        self.collection.add_fo_component(self)

    def _init_elements(self, element_name, element_number, collection):
        name_coarse = f'fo_coarse_{element_name}_{element_number}'
        self.fo_line_coarse = FanOutLineCoarse(name_coarse, collection)
        name_fine = f'fo_line_{element_name}_{element_number}'
        self.fo_line_fine = FanOutLineFine(name_fine, collection)
        self.via = None
        self.via_etch = None

    def add_via(self):
        if self.via is not None:
            self.via.layer = self.layer
            self.via.layer_stage = 'via_fine'
            self.via.build()

            uc_via = self.add_component(self.via)
            center = self.fo_points.qda.elements[self.element_name]['positions'][self.element_number]
            uc_via.center = center
            uc_via.build()

    def add_via_etch(self):
        if self.via is not None:
            self.via_etch.layer = self.layer
            self.via_etch.layer_stage = 'via_etch'
            self.via_etch.build()

            uc_via_etch = self.add_component(self.via_etch)
            center = self.fo_points.qda.elements[self.element_name]['positions'][self.element_number]
            uc_via_etch.center = center
            uc_via_etch.build()

    def add_coarse_fo_line(self):

        # self.fo_line_coarse.name = name
        # self.fo_line_coarse.cell.name = name
        self.elements[self.fo_line_coarse.name] = {'vertices': [],
                                                   'positions': [],
                                                   'layer': self.fo_line_coarse.layer}
        # self.fo_line_coarse = self.fo_line_coarse(name)

        self.fo_line_coarse.fo_direction = self.fo_direction

        self.fo_line_coarse.layer = self.layer
        if self.via is not None:
            self.fo_line_coarse.layer_stage = 'via_coarse'
        else:
            self.fo_line_coarse.layer_stage = 'coarse'

        attributes = {'fo_direction': self.fo_direction,
                      'n_fanout': self.n_fanout,
                      'fo_points': self.fo_points}
        missing_attrs = [name for name,
                         value in attributes.items() if value is None]
        if missing_attrs:
            raise ValueError(
                f"Attributes {', '.join(missing_attrs)} cannot be None.")
        if isinstance(self.element, Ohmic):
            fo_poly = self.fo_points.fo_ohmics_polygons_coarse[self.fo_direction][self.n_fanout]
        else:
            fo_poly = self.fo_points.fo_polygons_coarse[self.fo_direction][self.n_fanout]
        self.fo_line_coarse.polygons = fo_poly

        # self.elements[name] = self.fo_line_coarse
        self.add_component(self.fo_line_coarse, build=True)

    def add_fine_fo_line(self):
        # name = f'fo_line_{self.element_name}_{self.element_number}_fine'
        # self.fo_line_fine.name = name
        # self.fo_line_fine.cell.name = name
        if self.fo_direction_start:
            contact = [[0.01*self.fo_direction_start[0],
                        0.01*self.fo_direction_start[1],
                        'start']]
            self.fo_line_fine.points_along_path = contact + \
                self.fo_line_fine.points_along_path

        self.elements[self.fo_line_fine.name] = {'vertices': [],
                                                 'positions': [],
                                                 'layer': self.fo_line_fine.layer}

        self.fo_line_fine.fo_direction = self.fo_direction

        self.fo_line_fine.layer = self.layer
        if self.via is not None:
            self.fo_line_fine.layer_stage = 'via_fine'
        else:
            self.fo_line_fine.layer_stage = 'fine'

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
        elif self.element_name not in self.fo_points.qda.elements:
            raise KeyError(f"'{self.element_name}' not found in QuantumDotArray elements. "
                           "Have you built the QuantumDotArray before?")
        else:
            self.fo_line_fine.fo_start = self.fo_points.qda.elements[
                self.element_name]['positions'][self.element_number]
            self.fo_line_fine.fo_start[0] += self.start_offset[0]
            self.fo_line_fine.fo_start[1] += self.start_offset[1]
            self.fo_line_fine.fo_end = self.fo_points.get_fo_overlap_points(
                self.n_fanout, self.fo_direction)

        # self.elements[name] = self.fo_line_fine
        self.add_component(self.fo_line_fine, build=True)

    def plot(self):
        fig, axs = plt.subplots(1, 2, figsize=(6, 3))

        super().plot(axs[0])
        super().plot(axs[1])

        axs[0].set_xlim(-1.1*self.fo_points.fo_stages[0][0],
                        1.1*self.fo_points.fo_stages[0][0])
        axs[0].set_ylim(-1.1*self.fo_points.fo_stages[0][1],
                        1.1*self.fo_points.fo_stages[0][1])

        axs[0].set_title('fine fanout')
        axs[1].set_title('coarse fanout')

        fig.suptitle(self.name)

        plt.tight_layout()

    def build(self):
        self.add_via_etch()
        self.add_via()
        self.add_coarse_fo_line()
        self.add_fine_fo_line()
        super().build()
