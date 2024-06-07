# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:08:32 2023

@author: vjohn
"""
# %% imports

import numpy as np
import matplotlib.pyplot as plt

from ..helpers.helpers import apply_sublattice

from abc import ABC

# %% definition


class PlotMixin(ABC):
    def plot(self, ax=None):
        # If no axis is provided, create one
        if ax is None:
            fig, ax = plt.subplots()

        # Collect unique layers
        layer_numbers = []
        for el in self.elements.values():
            layer = el['layer']
            layer_stage = el['layer_stage']
            layer_number = getattr(layer, layer_stage)
            layer_numbers.append(layer_number)
        unique_layers = list(set(layer_numbers))

        # Create a colormap with enough colors for each layer
        colors = plt.cm.jet(np.linspace(0, 1, len(unique_layers)))

        # Create a dictionary to map layer values to colors
        layer_to_color = dict(zip(unique_layers, colors))

        all_x_positions = []  # Collect all x positions to set the x-axis limits
        all_y_positions = []  # Collect all y positions to set the y-axis limits

        max_polygon_size = 0  # To store the maximum size of the polygon

        # Plot each element
        for _, el in self.elements.items():
            vertices = el['vertices']
            positions = el['positions']

            combined_vertices = apply_sublattice(positions, vertices)

            # Calculate the size of the polygon
            x_positions, y_positions = zip(*vertices)
            polygon_width = max(x_positions) - min(x_positions)
            polygon_height = max(y_positions) - min(y_positions)
            polygon_size = max(polygon_width, polygon_height)
            max_polygon_size = max(max_polygon_size, polygon_size)

            for i in range(0, len(combined_vertices), len(vertices)):
                polygon_vertices = combined_vertices[i:i+len(vertices)]
                layer_number = getattr(el['layer'], el['layer_stage'])
                polygon = plt.Polygon(
                    polygon_vertices, color=layer_to_color[layer_number])
                ax.add_patch(polygon)

                # Extract x and y coordinates for axes limits
                x_positions, y_positions = zip(*polygon_vertices)
                all_x_positions.extend(x_positions)
                all_y_positions.extend(y_positions)

        # Set the axes limits with the maximum size of the polygon as margin
        ax.set_xlim(min(all_x_positions) - max_polygon_size,
                    max(all_x_positions) + max_polygon_size)
        ax.set_ylim(min(all_y_positions) - max_polygon_size,
                    max(all_y_positions) + max_polygon_size)

        ax.set_xlabel(r'width ($\mathrm{\mu m}$)')
        ax.set_ylabel(r'height ($\mathrm{\mu m}$)')
        ax.set_title(self.name)
        ax.set_aspect('equal', adjustable='box')

        if ax is None:
            plt.show()
