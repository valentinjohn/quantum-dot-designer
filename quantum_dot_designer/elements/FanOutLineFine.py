# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from ..elements.FanOutLineBase import FanOutLineBase
import numpy as np
from ..helpers.helpers import get_polygons_from_path
from ..BaseCollection import BaseCollection

class FanOutLineFine(FanOutLineBase):
    def __init__(self, name, collection: BaseCollection):
        super().__init__(name, collection)
        self._layer_stage = 'fine'
        self.fo_width_start = 40e-3
        self.fo_start = None
        self.fo_end = None
        self.points_along_path = []

    def calculate_fine_fo(self):
        """
        Calculate a refined path or polygon based on given points.

        Args:
        - start (list): Starting point in the format [x, y, width].
        - end (list): Ending point in the format [x, y, width].
        - points_along_path (list): List of points along the path. Each point can optionally have width and reference type.
        - return_type (str): Determines what to return. 'polygon' returns the vertices of the polygon.
                              'path' returns the absolute path points with their widths.

        Returns:
        - list: Vertices of the polygon if return_type is 'polygon'. Path points with their widths if return_type is 'path'.
        """

        start = self.fo_start.copy()
        start.append(self.fo_width_start)
        end = self.fo_end[1]
        points_along_path = self.points_along_path.copy()
        if self.fo_direction in ['left', 'right']:
            points_along_path.append([0.95, 1, self.fo_end[0][2], 'dif'])
        else:
            points_along_path.append([1, 0.95, self.fo_end[0][2], 'dif'])
        points_along_path.append(self.fo_end[0])

        # points_along_path

        path_points = [start]  # Starting with the start point
        current_point = np.array(start[:2])
        prev_width = float(start[2])
        for point in points_along_path:
            x, y = point[0], point[1]
            if len(point) == 3:
                # If the third element is a float, it's the width
                if isinstance(point[2], (float, int)):
                    width = float(point[2])
                    reference = 'abs'
                # Otherwise, it's the reference
                else:
                    reference = point[2]
                    width = prev_width
            elif len(point) == 4:
                # Check which entry is the float and which is the string
                if isinstance(point[2], (float, int)):
                    width = float(point[2])
                    reference = point[3]
                else:
                    width = point[3]
                    reference = point[2]
            else:
                width = prev_width
                reference = 'abs'

            prev_width = width  # Update the previous width

            if reference in ['absolute', 'abs', 'a']:
                next_point = np.array([x, y])
            elif reference in ['relative_to_start', 'start', 's']:
                next_point = np.array(start[:2]) + [x, y]
            elif reference in ['relative_to_previous', 'prev', 'p']:
                next_point = current_point + [x, y]
            elif reference in ['relative_to_difference', 'dif', 'd']:
                delta = np.array(self.fo_end[0][:2]) - np.array(start[:2])
                next_point = np.array(start[:2]) + delta * np.array([x, y])
            else:
                raise ValueError(f"Unknown reference type: {reference}")

            path_points.append([next_point[0], next_point[1], width])
            current_point = next_point

        path_points.append(end)  # Add the end point to the path

        self.path = path_points
        poly_vertices = get_polygons_from_path(path_points)
        self.polygons = poly_vertices

        return poly_vertices

    def build(self):
        self.calculate_fine_fo()
        super().build()
