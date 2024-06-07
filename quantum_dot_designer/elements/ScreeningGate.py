# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""
from ..base import Element, Layer
from ..BaseCollection import BaseCollection
import numpy as np
from ..helpers.helpers import (create_segmented_path,
                                                  get_polygons_from_path,
                                                  get_end_path,
                                                  point_along_line,
                                                  orthogonal_unit_vector,
                                                  adjust_vector_direction)

import gdstk


class ScreeningGate(Element):
    """
    Represents a screening gate within the quantum_dot_designer system, used for controlling electrostatic potentials in quantum dot device designs.

    The ScreeningGate class extends the standard Element class, adding functionalities that allow it to manage and manipulate screening paths. These paths are crucial in the design of quantum dot devices, as they help control the electrostatic environments within these devices, a key factor in their operational behavior.

    This class maintains a collection of paths and vertices representing the screening effect and provides methods to create segmented paths and retrieve geometric representations from these paths. The build process involves constructing these paths and registering their details within the device's design specifications.

    Attributes:
        layer (Layer): The layer to which the screening gate belongs. It must be a valid Layer instance.
        _layer_stage (str): The stage of the layer associated with the screening gate, defaults to 'fine'.
        vertices (list): A collection of vertices defining the geometry of the screening paths.
        screen_paths (list): A collection of segmented paths representing the screening effect.
        collection (BaseCollection): The collection of elements and paths to which the screening gate belongs.
        _screen_path (bool): A flag indicating whether a screening path is active or defined.
        contact_vertices (list): Vertices defining the contact points of the screening paths.
        fo_contact_width (float): The width of the contact point along the screening path.
        fo_contact_direction (int): The direction indicator for the contact point along the screening path.

    Methods:
        screen(element_name, element_number, points, widths): Create a segmented path based on specified parameters and add it to the screen paths collection.
        build(): Constructs the geometric representation of the screening gate and its paths, and integrates them into the quantum dot device design.
    """

    def __init__(self, name, collection: BaseCollection):
        super().__init__(name, collection)
        self.layer = None
        self._layer_stage = 'fine'
        self.vertices = []
        self.screen_paths = []
        self.collection = collection
        self._screen_path = False
        self.contact_vertices = None
        self.fo_contact_width = 40e-3
        self.fo_contact_direction = 0

    def screen(self, element_name, element_number,
               points, widths):
        """
        Create a segmented path for screening, based on the parameters of an existing element.

        This method generates a new segmented path that represents a screening effect, based on the specified element's properties. It calculates the path's geometry, updates the screen paths collection, and stores the vertices for later use in the build process.

        Args:
            element_name (str): The name of the element based on which the screening path is created.
            element_number (int): The identifier number of the element.
            points (list): A list of points that define the segments of the path.
            widths (list): A list of widths corresponding to each segment of the path.
        """
        fo_line_name = f'fo_line_{element_name}_{element_number}'
        element_fo_path = self.collection.fo_elements[fo_line_name].path
        screen_path = create_segmented_path(element_fo_path,
                                            points,
                                            widths)
        self.screen_paths.append(screen_path)
        poly_vertices = get_polygons_from_path(screen_path)
        self.vertices.append(poly_vertices)

        self._screen_path = True

    def build(self):
        """
        Construct the screening gate's geometric representation and integrate it into the device design.

        This method processes the stored vertices and screen paths, creating a detailed geometric representation of the screening gate. It handles the construction of the actual paths, applies necessary geometric transformations, and registers the details within the device's design specifications.

        The build process ensures that all screening paths are accurately represented and that their geometric details are correctly integrated into the overall quantum dot device design.

        Raises:
            ValueError: If no valid Layer is assigned before building.
        """
        if self.layer is None or not isinstance(self.layer, Layer):
            raise ValueError(
                "A valid Layer must be assigned before building the element.")

        cell = gdstk.Cell(self.name)
        layer = getattr(self.layer, self._layer_stage)
        for vertices in self.vertices:
            screen = gdstk.Polygon(vertices,
                                   layer=layer)
            screen.fillet(self.fillet, tolerance=self.fillet_tolerance)
            cell.add(screen)
            self.elements[self.name]['vertices'] = screen.points
            self.elements[self.name]['positions'] = [[0, 0]]
            self.elements[self.name]['layer'] = self.layer
            self.elements[self.name]['layer_stage'] = self._layer_stage

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
