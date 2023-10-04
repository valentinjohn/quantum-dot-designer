# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 10:05:40 2023

@author: vjohn
"""

from QuantumDotDesigner import QuantumDotArray
from QuantumDotDesigner.helpers.helpers import (compute_fanout_positions,
                                                get_fo_lines,
                                                generate_polygon_for_fanout)


class FanoutPoints():
    def __init__(self, qda):
        if not isinstance(qda, QuantumDotArray):
            raise TypeError(
                f"Expected qda to be of type {QuantumDotArray}, but got {type(qda)} instead.")
        if not qda.built:
            qda.build()
        self.qda = qda
        self.fo_lines = {}
        self._all_directions = ['top', 'bottom', 'right', 'left']
        self.fo_polygons_coarse = None
        self.fo_stages = [(16, 16), (500, 530), (1200, 1200)]
        self.fo_widths = [1, 6, 25]
        self.fanout_counts = {'top': 14, 'bottom': 14, 'left': 13, 'right': 13}
        self.spacings = [2, 80, 200]
        self.fo_fine_coarse_overlap = 3
        self.bondpad_position = {'top':  1500, 'bottom': 1500,
                                 'left': 1500, 'right': 1500}
        self.bondpad_size = {'top': (110, 400), 'bottom': (110, 400),
                             'left': (400, 110), 'right': (400, 110)}
        self._n = 0
        self.fanout_positions = None

    def create_fo_polygons_coarse(self):
        polygons = {}
        fanout_positions = compute_fanout_positions(self.fo_stages,
                                                    self.fanout_counts,
                                                    self.spacings)
        self.fanout_positions = fanout_positions
        self.fo_lines = get_fo_lines(fanout_positions,
                                     self.fanout_counts, self.fo_stages)
        for direction in self._all_directions:
            polygons[direction] = [generate_polygon_for_fanout(direction,
                                                               self.fo_lines[direction][n_fo],
                                                               self.bondpad_position,
                                                               self.bondpad_size,
                                                               self.fo_widths,
                                                               self.fo_fine_coarse_overlap)
                                   for n_fo in range(self.fanout_counts[direction])]

        self.fo_polygons_coarse = polygons

    def get_fo_overlap_points(self, n_fanout, direction):
        multiplier = 1
        if direction in ['bottom', 'left']:
            multiplier = -1
        if direction in ['top', 'bottom']:
            overlap_start = multiplier * \
                (self.fo_stages[0][1] - self.fo_fine_coarse_overlap)
            overlap_end = multiplier*self.fo_stages[0][1]
            fanout_positions = self.fanout_positions['device'][direction][n_fanout]
            fo_end_width = self.fo_widths[0]

            fo_end_points = [[fanout_positions, overlap_start, fo_end_width],
                             [fanout_positions, overlap_end, fo_end_width]]
        else:
            overlap_start = multiplier * \
                (self.fo_stages[0][0] - self.fo_fine_coarse_overlap)
            overlap_end = multiplier*self.fo_stages[0][0]
            fanout_positions = - \
                self.fanout_positions['device'][direction][n_fanout]
            fo_end_width = self.fo_widths[0]

            fo_end_points = [[overlap_start, fanout_positions, fo_end_width],
                             [overlap_end, fanout_positions, fo_end_width]]
        return fo_end_points
