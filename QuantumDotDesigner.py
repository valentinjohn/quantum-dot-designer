# -*- coding: utf-8 -*-
"""
Created on Tue May  9 07:55:18 2023

This module defines classes and functions related to a quantum dot array.

Classes:
- QuantumDotArrayElements: Represents a collection of quantum dot array elements.
- UnitCell: Represents a unit cell in the quantum dot array.
- QuantumDotArray: Represents a quantum dot array.
- Element (abstract): Represents a quantum dot array element.
- Plunger: Represents a plunger element in the quantum dot array.
- Barrier: Represents a barrier element in the quantum dot array.
- ScreeningGate: Represents a screening gate element in the quantum dot array.
- Ohmic: Represents an ohmic element in the quantum dot array.
- Sensor: Represents a sensor in the quantum dot array.
- Sublattice: Represents a sublattice in the quantum dot array.

Functions:
- rot_mat: Computes the rotation matrix for a given angle.
- gen_poly: Generates a regular polygon with a given number of sides.

"""

# %% imports
from dataclasses import dataclass
# from __future__ import annotations
from typing import Optional
import math
import gdstk
import numpy as np
import copy
from abc import ABC, abstractmethod
from collections import defaultdict
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from helpers.helpers import *

# %% Base


class PlotMixin:
    def plot(self, ax=None, build=False):
        if not self._built:
            if build:
                self.build()
            else:
                raise Exception(
                    'Copmonent not built yet, and therefore cannot be plotted.')
        # If no axis is provided, create one
        if ax is None:
            fig, ax = plt.subplots()

        # Collect unique layers
        unique_layers = list({el['layer'] for el in self.elements.values()})

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
                polygon = plt.Polygon(
                    polygon_vertices, color=layer_to_color[el['layer']])
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


class Element(ABC, PlotMixin):
    def __init__(self, name):
        """
        Initialize an Element object.

        Args:
            name (str): Name of the element.
            layer (int): Description of layer. Default is None.
            rotate (float): Rotation of the element. Default is 0.0.
            fillet (float): Fillet value. Indicates how much the polygon is rounded.
        """

        self.name = name
        self.layer = None
        self.elements = {name: {'vertices': [],
                                'positions': [],
                                'layer': self.layer}}
        self.cell = None
        self.rotate = 0.0
        self.fillet = 0.0
        self.fillet_tolerance = 1e-3
        self._built = False

    # default attributes to skip for Element
    _skip_copy_attrs = {'name', 'cell', 'elements'}

    @property
    def built(self):
        return self._built

    # This setter is private and can only be used within the class
    def _set_built(self, new_value: bool):
        self._built = new_value

    def copy(self, copy_name):
        # Use the same class as the current instance
        copied_obj = self.__class__(copy_name)

        for attr, value in self.__dict__.items():
            if attr not in self._skip_copy_attrs:
                setattr(copied_obj, attr, value)

        return copied_obj

    @abstractmethod
    def build(self):
        """
        Build the element.

        This method should be implemented in subclasses to build the specific element.
        """
        pass


class Sublattice(PlotMixin):
    def __init__(self, name):
        """
        Initialize a Sublattice object.

        Args:
            name (str): Name of the sublattice.
        """
        self.name = name
        self.elements = {}
        self.component = None
        self.rows = 1
        self.columns = 1
        self.spacing = (100, 100)
        self.center = (0, 0)
        self.cell = None
        self.xlim = None
        self.ylim = None
        self._width = (self.columns-1) * self.spacing[0]
        self._height = (self.rows-1) * self.spacing[1]
        self._built = False

    @property
    def built(self):
        return self._built

    # This setter is private and can only be used within the class
    def _set_built(self, new_value: bool):
        self._built = new_value

    def _update_width(self):
        """
        Update the width of the sublattice based on the number of columns and
        spacing.
        """
        self._width = (self.columns-1) * self.spacing[0]

    def _update_height(self):
        """
        Update the height of the sublattice based on the number of rows and
        spacing.
        """
        self._height = (self.rows-1) * self.spacing[1]

    def _get_lim(self, axis=0):
        """
        Get the minimum and maximum limits along the specified axis of the
        sublattice cell.

        Args:
            axis (int): Axis along which to find the limits (0 for x-axis,
                                                             1 for y-axis).

        Returns:
            tuple: Minimum and maximum limits as a tuple (min, max).
        """
        cell_max = 0
        cell_min = 0
        for poly in self.cell.get_polygons():
            poly_max = poly.points[:, axis].max()
            poly_min = poly.points[:, axis].min()
            cell_max = max(cell_max, poly_max)
            cell_min = min(cell_min, poly_min)
        if axis == 0:
            self.xlim = (cell_min, cell_max)
        elif axis == 1:
            self.ylim = (cell_min, cell_max)
        else:
            raise ValueError(
                f"axis has to be 0 or 1, but is '{axis}'.")

    def build(self):
        """
        Build the sublattice cell by adding the element references to the
        sublattice cell.
        """
        if not self.component.built:
            self.component.build()

        self._update_height()
        self._update_width()
        cell = gdstk.Cell(self.name)
        cell.add(gdstk.Reference(self.component.cell,
                                 (self.center[0] - self._width/2,
                                  self.center[1] - self._height/2),
                                 columns=self.columns,
                                 rows=self.rows,
                                 spacing=(self.spacing[0], self.spacing[1])
                                 ))
        self.cell = cell
        self._set_built(True)
        self._get_lim(axis=0)
        self._get_lim(axis=1)
        # only update elements attribute if it exists
        if not isinstance(self.component, FanOutLineBase):
            self.elements = update_positions(self.component.elements,
                                             self.rows, self.columns,
                                             self.spacing[0], self.spacing[1],
                                             self.center)
        else:
            pass


class UnitCell(PlotMixin):
    def __init__(self, name):
        """
        Initialize a UnitCell object.

        Args:
            parent_instance: The QuantumDotArray instance that this UnitCell belongs to.
            name (str): Name of the unit cell (default: 'unit_cell').
        """

        self.name = name
        self.elements = {}
        self.components = {}
        self.cell = gdstk.Cell(name)
        self._xlim = None
        self._ylim = None
        self._n = 0
        self._built = False

    @property
    def built(self):
        return self._built

    # This setter is private and can only be used within the class
    def _set_built(self, new_value: bool):
        self._built = new_value

    def add_component(self, component=None, build=False):
        """
        Add a sublattice to the unit cell.

        Args:
            name (str): Name of the sublattice.
            component: You can assign otionally a component in the argument
            build: Adds and builds at the same time

        Returns:
            Sublattice: The created Sublattice object.
        """
        name = f'{self.name}_sublattice_{self._n}'
        self._n = self._n + 1
        sublattice = Sublattice(name)
        if component is not None:
            if component._built == False:
                component.build()
            sublattice.component = component
            if build:
                sublattice.build()
            else:
                pass
        else:
            pass
        self.components[name] = sublattice
        return sublattice

    def _get_lim(self, axis=0):
        cell_max = 0
        cell_min = 0
        for poly in self.cell.get_polygons():
            poly_max = poly.points[:, axis].max()
            poly_min = poly.points[:, axis].min()
            cell_max = max(cell_max, poly_max)
            cell_min = min(cell_min, poly_min)
        if axis:
            self._ylim = (cell_min, cell_max)
        else:
            self._xlim = (cell_min, cell_max)

    @property
    def xlim(self):
        if not self.built:
            self.build()
        return self._xlim

    @property
    def ylim(self):
        if not self.built:
            self.build()
        return self._ylim

    def build(self):
        """
        Build the unit cell by adding sublattices to the unit cell.

        This method adds the sublattices to the unit cell's cell object.
        """

        elements = {}
        for cell in self.components.values():
            if not cell._built:
                cell.build()
            self.cell.add(gdstk.Reference(cell.cell))
            elements = merge_device_positions(elements, cell.elements)
        self.elements = elements
        self.cell.flatten()
        self._get_lim(axis=0)
        self._get_lim(axis=1)
        self._set_built(True)


# %% Quantum Dot Array classes


class QuantumDotArrayElements:
    """
    Represents a collection of quantum dot array elements.

    Attributes:
        elements: Dictionary of quantum dot array elements.
        sublattices: List of sublattices in the quantum dot array.
        unit_cells: Dictionary of unit cells in the quantum dot array.

    Methods:
        add_plunger: Add a plunger element to the collection.
        add_barrier: Add a barrier element to the collection.
        add_screening_gate: Add a screening gate element to the collection.
        add_ohmic: Add an ohmic element to the collection.
        add_sensor: Add a sensor to the collection.

    """

    def __init__(self):
        self.components = {}

    def add_plunger(self, name: str):
        """
        Add a plunger element to the collection.

        Args:
            name: Name of the plunger element.

        Returns:
            The created plunger element.

        """
        plunger = Plunger(name)
        self.components[name] = plunger
        return plunger

    # def add(self, element):
    #     """
    #     Add a element to the collection.

    #     Args:
    #         name: Name of the element.

    #     Returns:
    #         The created plunger element.

    #     """
    #     self.components[name] = element

    def add_barrier(self, name: str):
        """
        Add a barrier element to the collection.

        Args:
            name: Name of the barrier element.

        Returns:
            The created barrier element.

        """
        barrier = Barrier(name)
        self.components[name] = barrier
        return barrier

    def add_screening_gate(self, name: str):
        """
        Add a screening gate element to the collection.

        Args:
            name: Name of the screening gate element.

        Returns:
            The created screening gate element.

        """
        screening_gate = ScreeningGate(name)
        self.components[name] = screening_gate
        return screening_gate

    def add_ohmic(self, name: str):
        """
        Add an ohmic element to the collection.

        Args:
            name: Name of the ohmic element.

        Returns:
            The created ohmic element.

        """
        ohmic = Ohmic(name)
        self.components[name] = ohmic
        return ohmic

    def add_sensor(self, name):
        """
        Add a sensor to the collection.

        Args:
            name: Name of the sensor.

        Returns:
            The created sensor.

        """
        sensor = Sensor(name, self)
        self.components[name] = sensor
        return sensor

    def add_clavier_gate(self, name):
        """
        Add a clavier to the collection.

        Args:
            name: Name of the clavier.

        Returns:
            The created clavier.

        """
        clavier_gate = ClavierGate(name)
        self.components[name] = clavier_gate
        return clavier_gate

    def add_clavier(self, name):
        """
        Add a clavier to the collection.

        Args:
            name: Name of the clavier.

        Returns:
            The created clavier.

        """
        clavier = Clavier(name, self)
        self.components[name] = clavier
        return clavier

    def add_fo_line_fine(self, name):
        fanout_line = FanOutLineFine(name)
        self.components[name] = fanout_line
        return fanout_line

    def add_fo_line_coarse(self, name):
        fanout_line = FanOutLineCoarse(name)
        self.components[name] = fanout_line
        return fanout_line

    def add_fo_line(self, element_name, element_number):
        fanout_line = FanOutLine(element_name, element_number, self)
        self.components[f'fo_{element_name}_{element_number}'] = fanout_line
        return fanout_line

    def add_copy(self, component, copy_name):
        # if not component.built:
        #     component.build()
        attributes = copy.copy(vars(component))
        attributes.pop('name')
        attributes.pop('cell')
        attributes.pop('elements')

        if 'components' in attributes:
            attributes.pop('components')
        # Check if 'qda_elements' is an attribute of the component
        if 'qda_elements' in attributes:
            qda_elements = attributes.pop('qda_elements')
            # Create new element using both 'copy_name' and 'qda_elements'
            new_element = type(component)(copy_name, qda_elements)
        else:
            new_element = type(component)(copy_name)
        new_element.__dict__.update(attributes)

        if isinstance(component, Sensor):
            new_element.cell = new_element.cell.copy(copy_name)
            new_element.plunger = component.plunger.copy(
                f'{copy_name}_plunger')
            new_element.barrier_source = component.barrier_source.copy(
                f'{copy_name}_barrier_source')
            new_element.barrier_drain = component.barrier_drain.copy(
                f'{copy_name}_barrier_drain')
            new_element.source = component.source.copy(
                f'{copy_name}_source')
            new_element.drain = component.drain.copy(
                f'{copy_name}_drain')
            new_element.barrier_sep = component.barrier_sep.copy(
                f'{copy_name}_barrier_seperation')

        if isinstance(component, Clavier):
            new_element.cell = new_element.cell.copy(copy_name)
            new_element.screen = component.screen.copy(
                f'{copy_name}_screen_clav')
            for element_name, clav_gate in component.clavier_gates.copy().items():
                new_element.clavier_gates[copy_name] = clav_gate.copy(
                    f'{copy_name}_{element_name}')

        self.components[copy_name] = new_element
        return new_element


class QuantumDotArray(PlotMixin):
    def __init__(self):
        """
        Initialize a QuantumDotArray object.

        Args:
            parent_instance: The parent instance.
            unitcell_instance (UnitCell): An optional UnitCell instance to use
            for the array.
        """
        self.name = 'Quantum_Dot_Array'
        self.spacing_qd = 200e-3
        self.spacing_qd_diag = 2**0.5 * self.spacing_qd
        self.elements = {}
        self.components = {}
        self.components_position = {}
        self.main_cell = gdstk.Cell('MAIN')
        self._n = 0
        self.chip_layout_path = "chip_layout.gds"
        self._built = False

    @property
    def built(self):
        return self._built

    def add_component(self):
        """
        Add a sublattice to the QuantumDotArray.

        Args:
            name (str): Name of the sublattice.

        Returns:
            Sublattice: The created Sublattice object.
        """
        name = f'MAIN_sublattice_{self._n}'
        self._n = self._n + 1
        sublattice = Sublattice(name)
        self.components[name] = sublattice
        return sublattice

    def add_chip_layout(self):
        layout = gdstk.read_rawcells(self.chip_layout_path)

        self.main_cell.add(gdstk.Reference(layout['TOP']))
        self.main_cell.flatten()

        return layout

    def build(self):
        """
        Build the QuantumDotArray.

        This method adds the sublattices and unit cells to the main cell object.
        """
        elements = {}
        for cell in self.components.values():
            if not cell.built:
                cell.build()
            elements = merge_device_positions(elements, cell.elements)
            if isinstance(cell, Sublattice):
                self.main_cell.add(gdstk.Reference(cell.cell))
            elif isinstance(cell, UnitCell):
                for c in cell.cells:
                    self.main_cell.add(gdstk.Reference(c.cell))
        self.main_cell.flatten()
        self.elements = elements
        self._built = True

    def save_as_gds(self, filename):
        """
        Save the QuantumDotArray as a GDS file.

        Args:
            filename (str): Name of the output GDS file.
        """
        lib = gdstk.Library()
        lib.add(self.main_cell, *self.main_cell.dependencies(True))
        lib.write_gds(filename)


# %% Elements


class Plunger(Element):
    def __init__(self, name):
        """
        Initialize an Plunger object.

        Args:
            name (str): Name of the element.
            layer (int): Description of layer. Default is None.
            rotate (float): Rotation of the plunger. Default is 0.0.
            fillet (float): Fillet value. Indicates how much the polygon is rounded.
            ... (other attributes)
        """
        super().__init__(name)
        self.layer = 21
        self.diameter = None
        self._asymx = 1
        self._asymy = 1 / self._asymx

    @property
    def asymx(self):
        return self._asymx

    @asymx.setter
    def asymx(self, value):
        if value != 0:  # to avoid division by zero
            self._asymx = value
            self._asymy = 1 / self._asymx

    @property
    def asymy(self):
        return self._asymy

    @asymy.setter
    def asymy(self, value):
        if value != 0:  # to avoid division by zero
            self._asymy = value
            self._asymx = 1 / self._asymy

    def build(self):
        """
        Build the plunger element.
        """
        pl_points = gen_poly(8)
        pl = gdstk.Polygon(pl_points, layer=self.layer)
        pl.scale(0.5 / np.cos(np.pi / 8) * self.diameter)
        pl.scale(sx=self._asymx, sy=self._asymy)
        # pl.translate(self.x, self.y)
        pl.fillet(self.fillet, tolerance=self.fillet_tolerance)
        pl.fillet(0.02, tolerance=1e-4)
        cell = gdstk.Cell(self.name)
        cell.add(pl)
        self.elements[self.name]['vertices'] = pl.points
        # self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.cell = cell
        self._set_built(True)


class Barrier(Element):
    def __init__(self, name):
        """
        Initialize a Barrier object.

        Args:
            name (str): Name of the barrier.
        """
        super().__init__(name)
        self.layer = 5
        self.width = None
        self.length = None

    def build(self):
        """
        Build the barrier element.
        """
        bar = gdstk.Polygon([(-self.length/2, -self.width/2),
                             (self.length/2, -self.width/2),
                             (self.length/2 + self.width/4, -self.width/8),
                             (self.length/2 + self.width/4, self.width/8),
                             (self.length/2, self.width/2),
                             (-self.length/2, self.width/2),
                             (-self.length/2 - self.width/2, self.width/8),
                             (-self.length/2 - self.width/2, -self.width/8)],
                            layer=self.layer)
        bar.rotate(self.rotate)
        # bar.translate(self.x, self.y)
        bar.fillet(self.fillet, tolerance=self.fillet_tolerance)
        cell = gdstk.Cell(self.name)
        cell.add(bar)
        self.elements[self.name]['vertices'] = bar.points
        # self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.cell = cell
        self._set_built(True)


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


class Ohmic(Element):
    def __init__(self, name):
        super().__init__(name)
        self.layer = 1
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
        ohmic = gdstk.Polygon(self.vertices, layer=self.layer)
        ohmic.fillet(self.fillet, tolerance=self.fillet_tolerance)
        ohmic.rotate(self.rotate)
        cell = gdstk.Cell(self.name)
        cell.add(ohmic)
        self.elements[self.name]['vertices'] = ohmic.points
        # self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.cell = cell
        self._set_built(True)


class ClavierGate(Element):
    def __init__(self, name):
        """
        Initialize a Clavier gate object.

        Args:
            name (str): Name of the clavier gate
        """
        super().__init__(name)
        self.layer = 25
        self.width = 100
        self.length = 100*4*10
        self.gate_width = 100
        self.gate_length = 300
        self.n_clav_rep = 10
        self.spacing = 100*4
        self.shift = 0
        self.x = 0
        self.y = 0
        self.fillet = 0  # 0.02
        self.fillet_tolerance = 1e-4
        self.rotation = 0

    def build(self):
        """
        Build the clavier gate element.
        """
        cl_points = generate_clavier_gates(self.width, self.length,
                                           self.gate_width,
                                           self.gate_length,
                                           self.n_clav_rep, self.spacing,
                                           self.shift, (-self.length/2, 0),
                                           self.rotation)
        cl = gdstk.Polygon(cl_points, layer=self.layer)

        cl.translate(self.x, self.y)
        cl.fillet(self.fillet, tolerance=self.fillet_tolerance)
        cell = gdstk.Cell(self.name)
        cell.add(cl)
        self.elements[self.name]['vertices'] = cl.points
        self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['layer'] = self.layer
        self.cell = cell
        self._set_built(True)


class Gate:
    def __init__(self, name):
        self.name = name
        self.points = None
        self.position = None
        self.layer = None


class FanOutLineBase(PlotMixin):
    def __init__(self, name):
        self.name = name
        self.element_name = None
        self.element_number = None
        self.fo_fine_coarse_overlap = None
        self.fo_fine_coarse_overlap_gap = 0.3
        self.layer = None
        self.polygons = None
        self.path = None
        self.fo_direction = None
        self.n_fanout = None
        self.cell = gdstk.Cell(self.name)
        self.fillet = 0
        self.fillet_tolerance = 1e-3
        self.elements = {self.name: {'vertices': [],
                                     'positions': [],
                                     'layer': self.layer}}
        self._built = False

    @property
    def built(self):
        return self._built

    # This setter is private and can only be used within the class
    def _set_built(self, new_value: bool):
        self._built = new_value

    def build(self):
        fo_line = gdstk.Polygon(self.polygons, layer=self.layer)
        fo_line.fillet(self.fillet, tolerance=self.fillet_tolerance)
        fo_line.fillet(0.02, tolerance=1e-4)

        self.elements[self.name]['vertices'] = fo_line.points
        self.elements[self.name]['positions'] = [0, 0]
        self.elements[self.name]['layer'] = self.layer

        self.cell.add(fo_line)
        self._set_built(True)


class FanOutLineFine(FanOutLineBase):
    def __init__(self, name):
        super().__init__(name)
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


class FanOutLineCoarse(FanOutLineBase):
    pass


# %% Components

class Sensor(UnitCell):
    def __init__(self, name: str, qda_elements: QuantumDotArrayElements):
        """
        Initialize a Sensor object.

        Args:
            name (str): Name of the sensor.
        """
        super().__init__(name)
        self.qda_elements = qda_elements
        self.plunger = qda_elements.add_plunger(f'{name}_plunger')
        self.barrier_source = qda_elements.add_barrier(
            f'{name}_barrier_source')
        self.barrier_drain = qda_elements.add_barrier(f'{name}_barrier_drain')
        self.source = qda_elements.add_ohmic(f'{name}_source')
        self.drain = qda_elements.add_ohmic(f'{name}_drain')
        self.barrier_sep = qda_elements.add_barrier(
            f'{name}_barrier_seperation')
        self.gap_ohmic_pl = 40
        self.gap_sep = 40
        self.source_pos = None
        self.source_pos_angle = None
        self.drain_pos = None
        self.drain_pos_angle = None
        self.sep_pos = None
        self.sep_pos_angle = None
        self.bar_drain_end = 'clockwise'
        self.bar_source_end = 'counter-clockwise'
        self.bar_sep_end = 'clockwise'
        self.source_position_offset = (0, 0)
        self.drain_position_offset = (0, 0)
        self.bar_sou_position_offset = (0, 0)
        self.bar_dra_position_offset = (0, 0)
        self.bar_sharp_source = (0, 0)
        self.bar_sharp_drain = (0, 0)
        self.fillet = (0, 1e-3)
        self.__feature_gap = None
        self.sd_position = None
        self.bar_position = None
        self.sep_position = None
        self.components_position = {}
        self._bar_angle_dict = {'top': 0, 'right': np.pi/2,
                                'bottom': np.pi, 'left': -np.pi/2,
                                'top-right': np.pi/4,
                                'bottom-right': 3/4*np.pi,
                                'bottom-left': -3/4*np.pi,
                                'top-left': -1/4*np.pi}

    def _calculate_positions(self):
        if self.source_pos:
            self.source_pos_angle = self._bar_angle_dict[self.source_pos]
        else:
            if self.source_pos_angle == None:
                raise Exception(
                    f'Either {self}.source_pos with "top", "bottom", ... or {self}.source_pos_angle with an angle in rad has to be defined, but both are None.')
        if self.drain_pos:
            self.drain_pos_angle = self._bar_angle_dict[self.drain_pos]
        else:
            if self.drain_pos_angle == None:
                raise Exception(
                    f'Either {self}.drain_pos with "top", "bottom", ... or {self}.drain_pos_angle with an angle in rad has to be defined, but both are None.')
        if self.sep_pos:
            self.sep_pos_angle = self._bar_angle_dict[self.sep_pos]
        else:
            if self.sep_pos_angle == None:
                raise Exception(
                    f'Either {self}.sep_pos with "top", "bottom", ... or {self}.sep_pos_angle with an angle in rad has to be defined, but both are None.')

        (i, j) = (np.sin(self.source_pos_angle), np.cos(self.source_pos_angle))
        (m, n) = (np.sin(self.drain_pos_angle), np.cos(self.drain_pos_angle))
        (u, v) = (np.sin(self.sep_pos_angle), np.cos(self.sep_pos_angle))

        plunger = self.plunger
        bar_source = self.barrier_source
        bar_drain = self.barrier_drain
        source = self.source
        drain = self.drain

        self.sd_position = (
            (i*(plunger._asymx*plunger.diameter/2 +
                self.gap_ohmic_pl) +
             self.source_position_offset[0],
             j*(plunger._asymy*plunger.diameter/2 +
                self.gap_ohmic_pl) +
             self.source_position_offset[1]),
            (m*(plunger._asymx*plunger.diameter/2 +
                self.gap_ohmic_pl) +
             self.drain_position_offset[0],
             n*(plunger._asymy*plunger.diameter/2 +
                self.gap_ohmic_pl) +
             self.drain_position_offset[1]))

        self.bar_position = (
            (i*(plunger._asymx*plunger.diameter/2 +
                bar_source.width/2-self.__feature_gap) +
             self.bar_sou_position_offset[0],
             j*(plunger._asymy*plunger.diameter/2 +
                bar_source.width/2-self.__feature_gap) +
             self.bar_sou_position_offset[1]),
            (m*(plunger._asymx*plunger.diameter/2 +
                bar_drain.width/2-self.__feature_gap) +
             self.bar_dra_position_offset[0],
             n*(plunger._asymy*plunger.diameter/2 +
                bar_drain.width/2 - self.__feature_gap) +
             self.bar_dra_position_offset[1]))

        self.sep_pos_angleition = (u*(plunger._asymx*plunger.diameter/2+self.gap_sep/2),
                                   v*(plunger._asymy*plunger.diameter/2+self.gap_sep/2))

    def _set_barrier_properties(self):
        source = self.source
        drain = self.drain
        bar_source = self.barrier_source
        bar_drain = self.barrier_drain
        bar_sep = self.barrier_sep

        if self.bar_drain_end == 'clockwise':
            drain.sensor_pos = 'top'
            drain.ohmic_pos = 'right'
        else:
            drain.sensor_pos = 'bottom'
            drain.ohmic_pos = 'right'

        if self.bar_source_end == 'clockwise':
            source.sensor_pos = 'top'
            source.ohmic_pos = 'right'
        else:
            source.sensor_pos = 'bottom'
            source.ohmic_pos = 'right'

        bar_drain_angle_offset = 0
        multiplier_bar_drain = 1
        if self.bar_drain_end == 'clockwise':
            bar_drain_angle_offset = np.pi
            multiplier_bar_drain = -1

        bar_source_angle_offset = 0
        multiplier_bar_source = 1
        if self.bar_source_end == 'clockwise':
            bar_source_angle_offset = np.pi
            multiplier_bar_source = -1

        bar_sep_angle_offset = 0
        if self.bar_sep_end == 'clockwise':
            bar_sep_angle_offset = np.pi

        # if isinstance(self.source_pos_angle, str):
        #     bar_source.rotate = self._bar_angle_dict[self.source_pos_angle]
        #     source.rotate = (self._bar_angle_dict[self.source_pos_angle] +
        #                      multiplier_bar_source*np.pi/2)
        # else:
        bar_source.rotate = -self.source_pos_angle+bar_source_angle_offset
        source.rotate = (-self.source_pos_angle+bar_source_angle_offset +
                         multiplier_bar_source*np.pi/2)
        # if isinstance(self.drain_pos_angle, str):
        #     bar_drain.rotate = -self._bar_angle_dict[self.drain_pos_angle]
        #     drain.rotate = (self._bar_angle_dict[self.drain_pos_angle] -
        #                     multiplier_bar_drain*np.pi/2)
        # else:
        bar_drain.rotate = -self.drain_pos_angle+bar_drain_angle_offset
        drain.rotate = (-self.drain_pos_angle+bar_drain_angle_offset +
                        multiplier_bar_drain*np.pi/2)
        # if isinstance(self.sep_pos_angle, str):
        #     bar_sep.rotate = self._bar_angle_dict[self.sep_pos_angle]
        # else:
        bar_sep.rotate = -self.sep_pos_angle+bar_sep_angle_offset

    def _build_and_add_elements(self):
        self.plunger.build()
        self.barrier_source.build()
        self.barrier_drain.build()
        self.source.build()
        self.drain.build()
        self.barrier_sep.build()

        self.add_component(self.plunger, build=True)

        sl_bs = self.add_component(self.barrier_source)
        sl_bs.center = self.bar_position[0]
        sl_bs.build()

        sl_bd = self.add_component(self.barrier_drain)
        sl_bd.center = self.bar_position[1]
        sl_bd.build()

        sl_bsep = self.add_component(self.barrier_sep)
        sl_bsep.center = self.sep_pos_angleition
        sl_bsep.build()

        sl_source = self.add_component(self.source)
        sl_source.center = self.sd_position[0]
        sl_source.build()

        sl_drain = self.add_component(self.drain)
        sl_drain.center = self.sd_position[1]
        sl_drain.build()

    def build_elements(self):
        self.__feature_gap = self.barrier_source.width - self.gap_ohmic_pl

        self._calculate_positions()

        self._set_barrier_properties()

        self.components_position = {
            self.plunger.name: (0, 0),
            self.barrier_source.name: (0, 0),
            self.barrier_drain.name: (0, 0),
            self.source.name: (0, 0),
            self.drain.name: (0, 0),
            self.barrier_sep.name: (0, 0)
        }

        self._build_and_add_elements()

        return self.cell

    def build(self):
        self.build_elements()
        super().build()


class Clavier(UnitCell):
    def __init__(self, name, qda_elements):
        super().__init__(name)
        self.qda_elements = qda_elements
        self.clavier_gates = {}
        self.screen = None
        self.clav_dot_size = 100e-3
        self.clav_gate_gap = 20e-3
        self.clav_width = 200e-3
        self._n_clav_gates = 4
        self.clav_gap = [0.450, 0.650]
        self.clav_layers = [25, 26]
        self.n_clav_rep = 8

        self.clav_gate_width = self.clav_dot_size
        self._clav_gate_length = list(np.array(self.clav_gap) +
                                      self.clav_dot_size)
        self._clav_length = ((self.clav_gate_width +
                             self.clav_gate_gap) *
                             self._n_clav_gates *
                             (self.n_clav_rep - 1) +
                             self.clav_gate_width)
        self._clav_gate_spacing = ((self.clav_length -
                                   self.clav_gate_width) /
                                   (self.n_clav_rep-1))

        self._screen_length = (self.clav_length +
                               (self._n_clav_gates-1) /
                               self._n_clav_gates *
                               self.clav_gate_spacing)
        self.screen_width = 100
        self.screen_gap = 0
        self.screen_layer = 3
        self._screen_position = (self.clav_dot_size +
                                 self.screen_width +
                                 2*self.screen_gap)

        self._spacing = ((self._n_clav_gates-1) *
                         self.clav_gate_width +
                         self._n_clav_gates * self.clav_gate_gap)
        self.x = 0
        self.y = 0
        self.fillet = 0.02
        self.fillet_tolerance = 1e-4
        self.rotation = 0
        self.mirror = False

    @property
    def screen_position(self):
        self._screen_position = (self.clav_dot_size +
                                 self.screen_width +
                                 2*self.screen_gap)
        return self._screen_position

    @property
    def clav_length(self):
        self._clav_length = ((self.clav_gate_width +
                             self.clav_gate_gap) *
                             self._n_clav_gates *
                             (self.n_clav_rep - 1) +
                             self.clav_gate_width)
        return self._clav_length

    @property
    def clav_gate_length(self):
        self._clav_gate_length = list(np.array(self.clav_gap) +
                                      self.clav_dot_size +
                                      self.screen_gap +
                                      self.screen_width)
        return self._clav_gate_length

    @property
    def clav_gate_spacing(self):
        self._clav_gate_spacing = ((self.clav_length -
                                   self.clav_gate_width) /
                                   (self.n_clav_rep-1))
        return self._clav_gate_spacing

    @property
    def screen_length(self):
        self._screen_length = (self.clav_length +
                               (self._n_clav_gates-1) /
                               self._n_clav_gates *
                               self.clav_gate_spacing)
        return self._screen_length

    @property
    def spacing(self):
        self._spacing = ((self._n_clav_gates-1) *
                         self.clav_gate_width +
                         self._n_clav_gates * self.clav_gate_gap)
        return self._spacing

    def _initialize_clavier_gates(self):
        self._sl_clavier_gates = {}
        self._clav_layers_order = self.clav_layers + self.clav_layers[::-1]

    def _build_even_clavier_gates(self, n):
        name = f'{self.name}_gate_{2*n}'
        x = (self.x + (n - (self._n_clav_gates - 1) / 2) *
             self.clav_gate_spacing / self._n_clav_gates)
        y = (self.y + self.clav_gate_length[n % 2] - self.clav_dot_size/2)

        if self.mirror:
            self.rotation = 180
            y = - y

        self.clavier_gates[name] = self.qda_elements.add_clavier_gate(name)
        self.clavier_gates[name].layer = self._clav_layers_order[2*n]
        self.clavier_gates[name].width = self.clav_width
        self.clavier_gates[name].length = self.clav_length
        self.clavier_gates[name].gate_width = self.clav_gate_width
        self.clavier_gates[name].gate_length = self.clav_gate_length[n % 2]
        self.clavier_gates[name].n_clav_rep = self.n_clav_rep
        self.clavier_gates[name].spacing = self.spacing
        self.clavier_gates[name].rotation = self.rotation
        self.clavier_gates[name].fillet = self.fillet
        self.clavier_gates[name].fillet_tolerance = self.fillet_tolerance
        self.clavier_gates[name].build()

        sl_name = f'sublattice_{self.name}_gate_{2*n}'
        self._sl_clavier_gates[sl_name] = self.add_component()
        self._sl_clavier_gates[sl_name].component = self.clavier_gates[name]
        self._sl_clavier_gates[sl_name].center = (x, y)
        self._sl_clavier_gates[sl_name].build()

    def _build_odd_clavier_gates(self, n):
        name_odd = f'{self.name}_gate_{2*n+1}'
        x = (self.x + (n+1/2)*self.clav_gate_spacing / self._n_clav_gates)
        y = (self.y + self.clav_gate_length[n % 2] - self.clav_dot_size/2)
        if not self.mirror:
            y = - y

        self.clavier_gates[name_odd] = self.qda_elements.add_clavier_gate(
            name_odd)
        self.clavier_gates[name_odd].layer = self._clav_layers_order[2*n]
        self.clavier_gates[name_odd].width = self.clav_width
        self.clavier_gates[name_odd].length = self.clav_length
        self.clavier_gates[name_odd].gate_width = self.clav_gate_width
        self.clavier_gates[name_odd].gate_length = self.clav_gate_length[n % 2]
        self.clavier_gates[name_odd].n_clav_rep = self.n_clav_rep
        self.clavier_gates[name_odd].spacing = self.spacing
        self.clavier_gates[name_odd].rotation = 180 + self.rotation
        self.clavier_gates[name_odd].fillet = self.fillet
        self.clavier_gates[name_odd].fillet_tolerance = self.fillet_tolerance
        self.clavier_gates[name_odd].build()

        sl_name_odd = f'sublattice_{self.name}_gate_{2*n+1}'
        self._sl_clavier_gates[sl_name_odd] = self.add_component()
        self._sl_clavier_gates[sl_name_odd].component = self.clavier_gates[name_odd]
        self._sl_clavier_gates[sl_name_odd].center = (x, y)
        self._sl_clavier_gates[sl_name_odd].build()

    def _build_screening_gate(self):
        self.screen = self.qda_elements.add_screening_gate(
            f'{self.name}_screen')
        self.screen.layer = self.screen_layer
        self.screen.vertices = [[(self.x-self.screen_length/2,
                                 self.y+self.screen_width/2),
                                (self.x+self.screen_length/2,
                                 self.y+self.screen_width/2),
                                (self.x+self.screen_length/2,
                                 self.y-self.screen_width/2),
                                (self.x-self.screen_length/2,
                                 self.y-self.screen_width/2)]]
        self.screen.build()

        sl_clav_screen = self.add_component()
        sl_clav_screen.component = self.screen
        sl_clav_screen.center = (0, 0)
        sl_clav_screen.rows = 2
        sl_clav_screen.spacing = (0, self.screen_position)
        sl_clav_screen.build()

    def build(self):
        self._initialize_clavier_gates()

        for n in range(int(self._n_clav_gates/2)):
            self._build_even_clavier_gates(n)
            self._build_odd_clavier_gates(n)

        self._build_screening_gate()

        super().build()

        return self.cell


class FanOutLine(UnitCell):
    def __init__(self, element_name: str, element_number: int,
                 qda_elements: QuantumDotArrayElements):
        super().__init__(f'fo_{element_name}_{element_number}')
        self.qda_elements = qda_elements
        self.fo = None
        name_coarse = f'fo_coarse_{element_name}_{element_number}'
        self.fo_line_coarse = qda_elements.add_fo_line_coarse(name_coarse)
        name_fine = f'fo_line_{element_name}_{element_number}'
        self.fo_line_fine = qda_elements.add_fo_line_fine(name_fine)
        self.element_name = element_name
        self.element_number = element_number
        # self.fo_fine_coarse_overlap = None
        # self.fo_fine_coarse_overlap_gap = 0.3
        self.layer = qda_elements.components[element_name].layer
        self.polygons = None
        self.fo_direction = None
        self.n_fanout = None
        self.cell = gdstk.Cell(self.name)
        self.fillet = 0
        self.fillet_tolerance = 1e-3
        self.elements = {}
        self.start_offset = (0, 0)
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

        if not any(attr is None for attr in [self.fo_direction, self.n_fanout, self.fo]):
            self.fo_line_coarse.polygons = self.fo.fo_polygons_coarse[
                self.fo_direction][self.n_fanout]

        # self.components[name] = self.fo_line_coarse
        self.add_component(self.fo_line_coarse, build=True)

    def add_fine_fo_line(self):
        # name = f'fo_line_{self.element_name}_{self.element_number}_fine'
        # self.fo_line_fine.name = name
        # self.fo_line_fine.cell.name = name
        self.elements[self.fo_line_fine.name] = {'vertices': [],
                                                 'positions': [],
                                                 'layer': self.fo_line_fine.layer}

        self.fo_line_fine.fo_direction = self.fo_direction
        self.fo_line_fine.layer = self.qda_elements.components[self.element_name].layer

        if not any(attr is None for attr in [self.fo_direction, self.n_fanout, self.fo]):
            self.fo_line_fine.fo_start = self.fo.qda.elements[
                self.element_name]['positions'][self.element_number]
            self.fo_line_fine.fo_start[0] = (self.fo_line_fine.fo_start[0] +
                                             self.start_offset[0])
            self.fo_line_fine.fo_start[1] = (self.fo_line_fine.fo_start[1] +
                                             self.start_offset[1])
            self.fo_line_fine.fo_end = self.fo.get_fo_overlap_points(self.n_fanout,
                                                                     self.fo_direction)

        # self.components[name] = self.fo_line_fine
        self.add_component(self.fo_line_fine, build=True)

    def build(self):
        self.add_coarse_fo_line()
        self.add_fine_fo_line()
        super().build()

# %% Fanout


class FanoutGenerator():
    def __init__(self, name, qda):
        if not qda.built:
            qda.build()
        self.name = name
        self.qda = qda
        self.elements = {}
        self.cell = gdstk.Cell(name)
        self.components = {}
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
        self._built = False

    @property
    def built(self):
        return self._built

    def create_fo_polygons_coarse(self):
        polygons = {}
        fanout_positions = compute_fanout_positions(self.fo_stages,
                                                    self.fanout_counts,
                                                    self.spacings)
        self.fanout_positions = fanout_positions
        fo_lines = get_fo_lines(fanout_positions,
                                self.fanout_counts, self.fo_stages)
        for direction in self._all_directions:
            polygons[direction] = [generate_polygon_for_fanout(direction,
                                                               fo_lines[direction][n_fo],
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

    def add_component(self, component=None, build=False):
        """
        Add a sublattice to the unit cell.

        Args:
            name (str): Name of the sublattice.
            component: You can assign optionally a component in the argument
            build: Adds and builds at the same time with a single component placed at origin.

        Returns:
            Sublattice: The created Sublattice object.
        """
        name = f'{self.name}_sublattice_{self._n}'
        self._n = self._n + 1
        sublattice = Sublattice(name)
        if component is not None:
            if not component._built:
                component.build()
            sublattice.component = component
            if build:
                sublattice.build()
            else:
                pass
        else:
            pass
        self.components[name] = sublattice
        return sublattice

    def build(self):
        elements = {}
        for cell in self.components.values():
            if not cell._built:
                cell.build()
            self.cell.add(gdstk.Reference(cell.component.cell))
            elements = merge_device_positions(
                elements, cell.component.elements)
        self.elements = elements
        self.cell.flatten()
        self._built = True
