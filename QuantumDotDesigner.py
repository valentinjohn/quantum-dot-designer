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

# %% definitions


def rot_mat(theta):
    """
    Compute the rotation matrix for a given angle.

    Args:
        theta: Angle of rotation in radians.

    Returns:
        Rotation matrix as a numpy array.

    """
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])


def gen_poly(n, sp=None):
    """
    Generate a regular polygon with a given number of sides.

    Args:
        n: Number of sides of the polygon.
        sp: Starting point of the polygon. Default is None.

    Returns:
        List of tuples representing the vertices of the polygon.

    """
    mat = rot_mat(2*np.pi/n)
    if sp is None:
        sp = [np.cos(np.pi/n), np.sin(np.pi/n)]
    poly = [sp]
    for k in range(n):
        poly.append(np.dot(mat, poly[-1]))
    return [tuple(co) for co in poly]


def rotate_point(point, angle, origin=(0, 0)):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    x, y = point
    ox, oy = origin

    qx = ox + math.cos(angle) * (x - ox) - math.sin(angle) * (y - oy)
    qy = oy + math.sin(angle) * (x - ox) + math.cos(angle) * (y - oy)

    return qx, qy


def generate_clavier_gates(width, length, gate_width, gate_length,
                           num_clavier, spacing, shift, position, rotation=0):
    vertices = []

    # Convert rotation to radians
    rotation = math.radians(rotation)

    # Calculate the total space occupied by the claviers and the gaps
    total_clavier_space = (num_clavier * gate_width +
                           (num_clavier - 1) * spacing)

    # Calculate the start point of the claviers (centered and shifted)
    start_clavier = (length - total_clavier_space) / 2 + shift

    # Start point
    vertices.append(rotate_point((0 + position[0], 0 + position[1]), rotation))

    # Bottom edge of the rectangle
    for i in range(num_clavier + 1):
        # Before clavier
        point = (start_clavier + i * (gate_width + spacing) -
                 spacing / 2 + position[0], 0 + position[1])
        vertices.append(rotate_point(point, rotation))

        # If there is a clavier here
        if i < num_clavier:
            # Move to the clavier
            point = (start_clavier + i * (gate_width + spacing) +
                     position[0], 0 + position[1])
            vertices.append(rotate_point(point, rotation))

            # Traverse the clavier
            point = (start_clavier + i * (gate_width + spacing) +
                     position[0], -gate_length + position[1])
            vertices.append(rotate_point(point, rotation))

            point = (start_clavier + i * (gate_width + spacing) +
                     gate_width + position[0],
                     - gate_length + position[1])
            vertices.append(rotate_point(point, rotation))

            point = (start_clavier + i * (gate_width + spacing) +
                     gate_width + position[0], 0 + position[1])
            vertices.append(rotate_point(point, rotation))

            # Return to the rectangle
            point = (start_clavier + i * (gate_width + spacing) +
                     gate_width + spacing / 2 + position[0],
                     0 + position[1])
            vertices.append(rotate_point(point, rotation))

    # Right end of the rectangle
    point = (length + position[0], 0 + position[1])
    vertices.append(rotate_point(point, rotation))

    # Top edge of the rectangle
    point = (length + position[0], width + position[1])
    vertices.append(rotate_point(point, rotation))

    point = (0 + position[0], width + position[1])
    vertices.append(rotate_point(point, rotation))

    # Closing the polygon
    point = (0 + position[0], 0 + position[1])
    vertices.append(rotate_point(point, rotation))

    return vertices


def create_lattice_positions(rows, columns, spacing_x, spacing_y, center):
    """
    Generate lattice positions based on input parameters.
    """
    x_center, y_center = center
    positions = []

    for i in range(rows):
        for j in range(columns):
            x = (j * spacing_x + x_center) - (columns - 1) * spacing_x / 2
            y = (i * spacing_y + y_center) - (rows - 1) * spacing_y / 2
            positions.append([x, y])
    return positions


def apply_sublattice(main_lattice, sub_lattice):
    """
    Apply the sub_lattice to each position in the main_lattice.
    """
    result = []

    for main_pos in main_lattice:
        for sub_pos in sub_lattice:
            result.append([main_pos[0] + sub_pos[0], main_pos[1] + sub_pos[1]])

    return result


def update_positions(elements, rows, columns, spacing_x, spacing_y, center):
    """
    Update the positions of elements based on the given lattice parameters.
    """
    # Generate lattice positions
    lattice_positions = create_lattice_positions(
        rows, columns, spacing_x, spacing_y, center)

    for key, value in elements.items():
        # Update the positions of each element in the dictionary
        value['positions'] = apply_sublattice(lattice_positions,
                                              value['positions'])

    return elements


def merge_dicts(dict1, dict2):
    """
    Merge two input dictionaries. If the names match but vertices or layers are different,
    raise an error. Otherwise, merge by combining the positions.
    """
    # Create a copy of the first dictionary to avoid modifying the original
    merged = dict1.copy()

    for key, value in dict2.items():
        # If the key (name) is already in the merged dictionary
        if key in merged:
            # Check if the vertices are the same
            if not np.array_equal(merged[key]['vertices'], value['vertices']):
                raise ValueError(
                    f"The vertices for key '{key}' are different between the two dictionaries.")

            # Check if the layers are the same
            if merged[key]['layer'] != value['layer']:
                raise ValueError(
                    f"The layers for key '{key}' are different between the two dictionaries.")

            # Merge the positions
            merged[key]['positions'] += value['positions']

            # Sort the positions based on y_max and x_min
            merged[key]['positions'].sort(key=lambda x: (-x[1], x[0]))

        else:
            # If the key (name) is not in the merged dictionary, just add the whole item
            merged[key] = value

    return merged


def compute_positions(length, n_lines, spacing):
    start = - (n_lines - 1) * spacing / 2
    return [start + i * spacing for i in range(n_lines)]


def compute_fanout_positions(rect_dims, fanout_counts, spacings):
    """
    Computes the positions of the fanout lines for each rectangle.

    Parameters:
    - rect_dims: A list of tuples, where each tuple contains the dimensions (length, width) of a rectangle.
    - fanout_counts: A dictionary with keys 'top', 'bottom', 'left', 'right' and values being the fanout counts.
    - spacings: A list of spacings for each rectangle.

    Returns:
    A dictionary containing the fanout positions for each side of each rectangle.
    """

    device_rect, intermediate_rect, bondpad_square = rect_dims
    device_spacing, intermediate_spacing, bondpad_spacing = spacings

    # Compute positions for each rectangle and side
    fanout_positions = {
        'device': {
            'top': compute_positions(device_rect[0], fanout_counts['top'], device_spacing),
            'bottom': compute_positions(device_rect[0], fanout_counts['bottom'], device_spacing),
            'left': compute_positions(device_rect[1], fanout_counts['left'], device_spacing),
            'right': compute_positions(device_rect[1], fanout_counts['right'], device_spacing)
        },
        'intermediate': {
            'top': compute_positions(intermediate_rect[0], fanout_counts['top'], intermediate_spacing),
            'bottom': compute_positions(intermediate_rect[0], fanout_counts['bottom'], intermediate_spacing),
            'left': compute_positions(intermediate_rect[1], fanout_counts['left'], intermediate_spacing),
            'right': compute_positions(intermediate_rect[1], fanout_counts['right'], intermediate_spacing)
        },
        'bondpad': {
            'top': compute_positions(bondpad_square[0], fanout_counts['top'], bondpad_spacing),
            'bottom': compute_positions(bondpad_square[0], fanout_counts['bottom'], bondpad_spacing),
            'left': compute_positions(bondpad_square[1], fanout_counts['left'], bondpad_spacing),
            'right': compute_positions(bondpad_square[1], fanout_counts['right'], bondpad_spacing)
        }
    }

    return fanout_positions


def get_fo_lines(fanout_positions, fanout_counts, rect_dims):
    fo_lines = {'top': [],
                'bottom': [],
                'right': [],
                'left': []}
    for n_fo in range(fanout_counts['top']):
        fo_line = [[fanout_positions['device']['top'][n_fo],
                    rect_dims[0][1]],
                   [fanout_positions['intermediate']['top'][n_fo],
                    rect_dims[1][1]],
                   [fanout_positions['bondpad']['top'][n_fo],
                    rect_dims[2][1]]
                   ]
        fo_lines['top'].append(fo_line)

    for n_fo in range(fanout_counts['bottom']):
        fo_line = [[fanout_positions['device']['bottom'][n_fo],
                    -rect_dims[0][1]],
                   [fanout_positions['intermediate']['bottom'][n_fo],
                    -rect_dims[1][1]],
                   [fanout_positions['bondpad']['bottom'][n_fo],
                    -rect_dims[2][1]]
                   ]
        fo_lines['bottom'].append(fo_line)

    for n_fo in range(fanout_counts['right']):
        fo_line = [[rect_dims[0][0],
                    -fanout_positions['device']['right'][n_fo]],
                   [rect_dims[1][0],
                    -fanout_positions['intermediate']['right'][n_fo]],
                   [rect_dims[2][0],
                    -fanout_positions['bondpad']['right'][n_fo]]
                   ]
        fo_lines['right'].append(fo_line)

    for n_fo in range(fanout_counts['left']):
        fo_line = [[-rect_dims[0][0],
                    -fanout_positions['device']['left'][n_fo]],
                   [-rect_dims[1][0],
                    -fanout_positions['intermediate']['left'][n_fo]],
                   [-rect_dims[2][0],
                    -fanout_positions['bondpad']['left'][n_fo]]
                   ]
        fo_lines['left'].append(fo_line)
    return fo_lines


def calculate_vertex(point, offset):
    """Helper function to add an offset to a given point and return the new point."""
    return list(np.array(point) + np.array(offset))


def generate_polygon_for_fanout(direction, fo_line, bondpad_position, bondpad_size, fo_widths, fo_fine_coarse_overlap):
    """Helper function to generate a polygon for a single fanout based on direction."""

    # Adjustments based on direction
    if direction in ['top', 'bottom']:
        coord_index = 0
        multiplier = 1 if direction == 'top' else -1
    else:  # 'left' or 'right'
        coord_index = 1
        multiplier = 1 if direction == 'right' else -1

    # Base points and offsets
    if direction in ['top', 'bottom']:
        base_point = [fo_line[2][0],
                      multiplier*abs(bondpad_position[direction])]
    else:
        base_point = [multiplier*abs(bondpad_position[direction]),
                      fo_line[2][1]]

    bondpad_offset = [multiplier * bondpad_size[direction]
                      [0]/2, multiplier * bondpad_size[direction][1]/2]
    fo_width_offset = [0, 0]
    fo_width_offset[coord_index] = multiplier * fo_widths[2]/2

    # Calculating bondpad corners
    if direction in ['top', 'bottom']:
        top_left = calculate_vertex(
            base_point, [-bondpad_offset[0], bondpad_offset[1]])
        top_right = calculate_vertex(
            base_point, [bondpad_offset[0], bondpad_offset[1]])
        bottom_left = calculate_vertex(
            base_point, [-bondpad_offset[0], -bondpad_offset[1]])
        bottom_right = calculate_vertex(
            base_point, [bondpad_offset[0], -bondpad_offset[1]])
        center_bottom_left = calculate_vertex(
            base_point, [-fo_widths[2]/2, -bondpad_offset[1]])
        center_bottom_right = calculate_vertex(
            base_point, [fo_widths[2]/2, -bondpad_offset[1]])

        bondpad = [center_bottom_right, bottom_right, top_right,
                   top_left, bottom_left, center_bottom_left]
        dir_unit = np.array([1, 0])
        dir_unit_orth = np.array([0, 1])
    else:
        top_left = calculate_vertex(
            base_point, [-bondpad_offset[0], bondpad_offset[1]])
        top_right = calculate_vertex(
            base_point, [bondpad_offset[0], bondpad_offset[1]])
        bottom_left = calculate_vertex(
            base_point, [-bondpad_offset[0], -bondpad_offset[1]])
        bottom_right = calculate_vertex(
            base_point, [bondpad_offset[0], -bondpad_offset[1]])
        center_bottom_left = calculate_vertex(
            base_point, [-bondpad_offset[0], -fo_widths[2]/2])
        center_top_left = calculate_vertex(
            base_point, [-bondpad_offset[0], fo_widths[2]/2])

        bondpad = [center_top_left, top_left, top_right,
                   bottom_right, bottom_left, center_bottom_left]
        dir_unit = np.array([0, 1])
        dir_unit_orth = np.array([1, 0])

    bondpad_to_fine = [calculate_vertex(fo_line[2], -dir_unit*fo_widths[2]/2),
                       calculate_vertex(fo_line[1], -dir_unit*fo_widths[1]/2),
                       calculate_vertex(fo_line[0], -dir_unit*fo_widths[0]/2)]

    overlap_fine_coarse = [calculate_vertex(fo_line[0],
                                            -dir_unit*fo_widths[0]/2 -
                                            dir_unit_orth*fo_fine_coarse_overlap * multiplier),
                           calculate_vertex(fo_line[0],
                                            dir_unit*fo_widths[0]/2 -
                                            dir_unit_orth*fo_fine_coarse_overlap * multiplier)]

    fine_to_bondpad = [calculate_vertex(fo_line[0], dir_unit*fo_widths[0]/2),
                       calculate_vertex(fo_line[1], dir_unit*fo_widths[1]/2),
                       calculate_vertex(fo_line[2], dir_unit*fo_widths[2]/2)]

    fo_polys = []
    fo_polys.extend(bondpad_to_fine)
    fo_polys.extend(overlap_fine_coarse)
    fo_polys.extend(fine_to_bondpad)
    fo_polys.extend(bondpad)

    return fo_polys
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

    def add_plunger(self, name):
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

    def add_barrier(self, name):
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

    def add_screening_gate(self, name):
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

    def add_ohmic(self, name):
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

    def add_copy(self, component, copy_name):
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
        self.components[copy_name] = new_element
        return new_element

    # def add_copy(self, component, copy_name):
    #     attributes = copy.copy(vars(component))
    #     attributes.pop('name')
    #     attributes.pop('cell')
    #     attributes.pop('elements')
    #     if 'components' in attributes:
    #         attributes.pop('components')
    #     # Check if 'qda_elements' is an attribute of the component
    #     if 'qda_elements' in attributes:
    #         qda_elements = attributes.pop('qda_elements')
    #         # Create new element using both 'copy_name' and 'qda_elements'
    #         new_element = type(component)(copy_name, qda_elements)
    #     else:
    #         new_element = type(component)(copy_name)
    #     new_element.__dict__.update(attributes)
    #     if isinstance(component, Sensor):
    #         new_element.plunger.name = f'{copy_name}_plunger'
    #         new_element.barrier_source.name = f'{copy_name}_barrier_source'
    #         new_element.barrier_drain.name = '{copy_name}_barrier_drain'
    #         new_element.sourcename = f'{copy_name}_source'
    #         new_element.drainname = f'{copy_name}_barrier_drain'
    #         new_element.barrier_sepname = f'{copy_name}_barrier_seperation'
    #     self.components[copy_name] = new_element
    #     return new_element

    def add_fo_line(self, element_name, n_element):
        name = f'fo_line_{element_name}_{n_element}'
        fo_line = FanOutLine(name)

        fo_line.element_name = element_name
        fo_line.n_element = n_element

        fo_line.layer_fine = self.components[element_name].layer
        fo_line.layer_coarse = fo_line.layer_fine + 20
        # fo_line.fog = self
        self.components[name] = fo_line

        return fo_line


class UnitCell:
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
        self.xlim = None
        self.ylim = None
        self._n = 0

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
            self.ylim = (cell_min, cell_max)
        else:
            self.xlim = (cell_min, cell_max)

    def build(self):
        """
        Build the unit cell by adding sublattices to the unit cell.

        This method adds the sublattices to the unit cell's cell object.
        """
        elements = {}
        for cell in self.components.values():
            self.cell.add(gdstk.Reference(cell.cell))
            elements = merge_dicts(elements, cell.elements)
        self.elements = elements
        self.cell.flatten()
        self._get_lim(axis=0)
        self._get_lim(axis=1)


class QuantumDotArray:
    def __init__(self):
        """
        Initialize a QuantumDotArray object.

        Args:
            parent_instance: The parent instance.
            unitcell_instance (UnitCell): An optional UnitCell instance to use
            for the array.
        """
        self.spacing_qd = 200e-3
        self.spacing_qd_diag = 2**0.5 * self.spacing_qd
        self.elements = {}
        self.components = {}
        self.components_position = {}
        self.main_cell = gdstk.Cell('MAIN')
        self._n = 0

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

    def build(self):
        """
        Build the QuantumDotArray.

        This method adds the sublattices and unit cells to the main cell object.
        """
        elements = {}
        for cell in self.components.values():
            elements = merge_dicts(elements, cell.elements)
            if isinstance(cell, Sublattice):
                self.main_cell.add(gdstk.Reference(cell.cell))
            elif isinstance(cell, UnitCell):
                for c in cell.cells:
                    self.main_cell.add(gdstk.Reference(c.cell))
        self.main_cell.flatten()
        self.elements = elements

    def save_as_gds(self, filename):
        """
        Save the QuantumDotArray as a GDS file.

        Args:
            filename (str): Name of the output GDS file.
        """
        lib = gdstk.Library()
        lib.add(self.main_cell)
        lib.write_gds(filename)

# %% Elements


class Element(ABC):
    def __init__(self, name):
        """
        Initialize an Element object.

        Args:
            name (str): Name of the element.
        """
        self.name = name
        self.layer = None
        self.elements = {name: {'vertices': [],
                                'positions': [],
                                'layer': self.layer}}
        self.cell = None
        self.x = 0
        self.y = 0
        self.rotate = 0.0
        self.fillet = 0
        self.fillet_tolerance = 1e-3

    @abstractmethod
    def build(self):
        """
        Build the element.

        This method should be implemented in subclasses to build the specific element.
        """
        pass


class Plunger(Element):
    def __init__(self, name):
        """
        Initialize a Plunger object.

        Args:
            name (str): Name of the plunger.
        """
        super().__init__(name)
        self.layer = 21
        self.diameter = None
        self.asym = 1.0
        self._asymx = self.asym
        self._asymy = 1 / self.asym

    def _update_asym(self, asym):
        """
        Update the asymmetry of the plunger.

        Args:
            asym (float): Asymmetry value.
        """
        self._asymx = asym
        self._asymy = 1 / asym

    def build(self):
        """
        Build the plunger element.
        """
        self._update_asym(self.asym)
        pl_points = gen_poly(8)
        pl = gdstk.Polygon(pl_points, layer=self.layer)
        pl.scale(0.5 / np.cos(np.pi / 8) * self.diameter)
        pl.scale(sx=self._asymx, sy=self._asymy)
        pl.translate(self.x, self.y)
        pl.fillet(self.fillet, tolerance=self.fillet_tolerance)
        pl.fillet(0.02, tolerance=1e-4)
        cell = gdstk.Cell(self.name)
        cell.add(pl)
        self.elements[self.name]['vertices'] = pl.points
        self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['layer'] = self.layer
        self.cell = cell


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
        bar.translate(self.x, self.y)
        bar.fillet(self.fillet, tolerance=self.fillet_tolerance)
        cell = gdstk.Cell(self.name)
        cell.add(bar)
        self.elements[self.name]['vertices'] = bar.points
        self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['layer'] = self.layer
        self.cell = cell


class ScreeningGate(Element):
    def __init__(self, name):
        super().__init__(name)
        self.layer = 5
        self.vertices = []

    def build(self):
        screen = gdstk.Polygon(self.vertices,
                               layer=self.layer)
        screen.translate(self.x, self.y)
        screen.fillet(self.fillet, tolerance=self.fillet_tolerance)
        cell = gdstk.Cell(self.name)
        cell.add(screen)
        self.elements[self.name]['vertices'] = screen.points
        self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['layer'] = self.layer
        self.cell = cell


class Ohmic(Element):
    def __init__(self, name):
        super().__init__(name)
        self.width = 0.0
        self.height = 0.0
        self.shape = None

    def build(self):
        print('Build method for Ohmic not implemeneted yet.')


class Sensor(UnitCell):
    def __init__(self, name: str, qda_elements: QuantumDotArrayElements):
        """
        Initialize a Sensor object.

        Args:
            name (str): Name of the sensor.
        """
        super().__init__(name)
        self.qda_elements = qda_elements
        # self.cell = None
        # self.elements = {name: {'vertices': [],
        #                         'positions': [],
        #                         'layer': self.layer}}
        self.x = 0
        self.y = 0
        self.plunger = qda_elements.add_plunger(f'{name}_plunger')
        self.barrier_source = qda_elements.add_barrier(
            f'{name}_barrier_source')
        self.barrier_drain = qda_elements.add_barrier(f'{name}_barrier_drain')
        self.source = Ohmic(f'{name}_source')
        self.drain = qda_elements.add_ohmic(f'{name}_barrier_drain')
        self.barrier_sep = qda_elements.add_barrier(
            f'{name}_barrier_seperation')
        self.gap_ohmic_pl = 40
        self.gap_sep = 40
        self.source_pos = 'left'
        self.drain_pos = 'right'
        self.sep_pos = 'bottom'
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

        # self.components = {self.plunger.name: self.plunger,
        #                    self.barrier_source.name: self.barrier_source,
        #                    self.barrier_drain.name: self.barrier_drain,
        #                    self.source.name: self.source,
        #                    self.drain.name: self.drain,
        #                    self.barrier_sep.name: self.barrier_sep
        #                    }

        self.components_position = {}

    def build_elements(self):
        """
        Build the sensor element.

        Returns:
            gdstk.Cell: The built sensor cell.
        """
        # cell = gdstk.Cell(self.name)
        plunger = self.plunger
        bar_source = self.barrier_source
        bar_drain = self.barrier_drain
        source = self.source
        drain = self.drain
        bar_sep = self.barrier_sep
        self.__feature_gap = self.barrier_source.width - self.gap_ohmic_pl

        orientation_dict = {'top': (0, 1), 'right': (1, 0),
                            'bottom': (0, -1), 'left': (-1, 0),
                            'top-right': (2**0.5/2, 2**0.5/2),
                            'bottom-right': (2**0.5/2, -2**0.5/2),
                            'bottom-left': (-2**0.5/2, -2**0.5/2),
                            'top-left': (-2**0.5/2, 2**0.5/2)}

        bar_angle_dict = {'top': np.pi, 'right': -np.pi/2,
                          'bottom': np.pi, 'left': +np.pi/2,
                          'top-right': -np.pi/4,
                          'bottom-right': np.pi/4,
                          'bottom-left': 3*np.pi/4,
                          'top-left': -3*np.pi/4}

        (i, j) = orientation_dict[self.source_pos]
        (m, n) = orientation_dict[self.drain_pos]
        (u, v) = orientation_dict[self.sep_pos]

        sd_position = ((i*(plunger._asymx*plunger.diameter/2 +
                           bar_source.width+source.width/2) +
                        self.source_position_offset[0],
                        j*(plunger._asymy*plunger.diameter/2 +
                           bar_source.width+source.width/2) +
                        self.source_position_offset[1]),
                       (m*(plunger._asymx*plunger.diameter/2 +
                           bar_source.width+source.width/2) +
                       self.drain_position_offset[0],
                       n*(plunger._asymy*plunger.diameter/2 +
                          bar_source.width+source.width/2) +
                       self.drain_position_offset[1]))

        bar_position = ((i*(plunger._asymx*plunger.diameter/2 +
                            bar_source.width/2-self.__feature_gap) +
                         self.bar_sou_position_offset[0],
                         j*(plunger._asymy*plunger.diameter/2 +
                            bar_source.width/2-self.__feature_gap) +
                         self.bar_sou_position_offset[1]),
                        (m*(plunger._asymx*plunger.diameter/2 +
                            bar_source.width/2-self.__feature_gap) +
                        self.bar_dra_position_offset[0],
                        n*plunger._asymy*(plunger.diameter/2 +
                                          bar_source.width/2 -
                                          self.__feature_gap) +
                        self.bar_dra_position_offset[1]))

        sep_position = (u*(plunger._asymx*plunger.diameter/2+self.gap_sep/2),
                        v*(plunger._asymy*plunger.diameter/2+self.gap_sep/2))

        self.sd_position = sd_position
        self.bar_position = bar_position

        bar_source.rotate = bar_angle_dict[self.source_pos]
        bar_source.x = bar_position[0][0]
        bar_source.y = bar_position[0][1]

        bar_drain.rotate = -bar_angle_dict[self.drain_pos]
        bar_drain.x = bar_position[1][0]
        bar_drain.y = bar_position[1][1]

        bar_sep.rotate = bar_angle_dict[self.sep_pos]
        bar_sep.x = sep_position[0]
        bar_sep.y = sep_position[1]

        self.components_position = {self.plunger.name:
                                    (self.plunger.x, self.plunger.y),
                                    self.barrier_source.name:
                                    (self.barrier_source.x,
                                     self.barrier_source.y),
                                    self.barrier_drain.name:
                                    (self.barrier_drain.x,
                                     self.barrier_drain.y),
                                    self.source.name:
                                    (self.source.x, self.source.y),
                                    self.drain.name:
                                    (self.drain.x, self.drain.y),
                                    self.barrier_sep.name:
                                    (self.barrier_sep.x, self.barrier_sep.y),
                                    }

        self.plunger.build()
        self.barrier_source.build()
        self.barrier_drain.build()
        # self.source.build()
        # self.drain.build()
        self.barrier_sep.build()

        self.add_component(self.plunger, build=True)
        self.add_component(self.barrier_source, build=True)
        self.add_component(self.barrier_drain, build=True)
        self.add_component(self.barrier_sep, build=True)

        # self.build()

        # cell.add(gdstk.Reference(plunger.cell))
        # cell.add(gdstk.Reference(bar_source.cell))
        # cell.add(gdstk.Reference(bar_drain.cell))
        # cell.add(gdstk.Reference(source))
        # cell.add(gdstk.Reference(drain))
        # cell.add(gdstk.Reference(bar_sep.cell))
        # self.cell = cell

        return self.cell


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
        self.fillet = 0.02
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
        self.cell = cell


class Clavier:
    def __init__(self, name, qda_elements):
        self.name = name
        self.qda_elements = qda_elements
        # self.qda = None
        self.cell = None
        self.clav_dot_size = 100
        self.clav_gate_gap = 20
        self.clav_width = 200
        self._n_clav_gates = 4
        self.clav_gap = [450, 650]
        self.clav_layers = [25, 26]
        self.n_clav_rep = 8

        self.clav_gate_width = self.clav_dot_size
        self.clav_gate_length = list(np.array(self.clav_gap) +
                                     self.clav_dot_size)
        self.clav_length = ((self.clav_gate_width +
                             self.clav_gate_gap) *
                            self._n_clav_gates *
                            (self.n_clav_rep - 1) +
                            self.clav_gate_width)
        self.clav_gate_spacing = ((self.clav_length -
                                   self.clav_gate_width) /
                                  (self.n_clav_rep-1))

        self.screen_length = (self.clav_length +
                              (self._n_clav_gates-1) /
                              self._n_clav_gates *
                              self.clav_gate_spacing)
        self.screen_width = 100
        self.screen_gap = 0
        self.screen_layer = 3

        self.spacing = ((self._n_clav_gates-1) *
                        self.clav_gate_width +
                        self._n_clav_gates * self.clav_gate_gap)
        self.x = 0
        self.y = 0
        self.fillet = 0.02
        self.fillet_tolerance = 1e-4

    def build(self):
        # if self.cell == None:
        #     raise Exception(f'QuantumDotDesigner.Clavier.qda is None. Assign' +
        #                     'QuantumDotDesigner.Clavier.qda a QuantumDotArray.')
        # else:
        #     pass

        cell = UnitCell()

        clavier_gates = {}
        sl_clavier_gates = {}

        clav_layers_order = self.clav_layers + self.clav_layers[::-1]

        for n in range(int(self._n_clav_gates/2)):
            clavier_gates[2 *
                          n] = self.qda_elements.add_clavier_gate(f'clavier_gates_{2*n}')
            clavier_gates[2*n].layer = clav_layers_order[2*n]
            clavier_gates[2*n].width = self.clav_width
            clavier_gates[2*n].length = self.clav_length
            clavier_gates[2*n].gate_width = self.clav_gate_width
            clavier_gates[2*n].gate_length = (self.clav_gate_length[n % 2]
                                              + n * self.clav_gate_width)
            clavier_gates[2*n].n_clav_rep = self.n_clav_rep
            clavier_gates[2*n].spacing = self.spacing
            clavier_gates[2*n].x = ((n - (self._n_clav_gates - 1) / 2) *
                                    self.clav_gate_spacing / self._n_clav_gates)
            clavier_gates[2*n].y = (self.y +
                                    self.clav_gate_length[n % 2]/2 +
                                    self.clav_gap[n % 2]/2 +
                                    n*self.clav_gate_width)
            clavier_gates[2*n].fillet = self.fillet
            clavier_gates[2*n].fillet_tolerance = self.fillet_tolerance
            clavier_gates[2*n].build()

            sl_clavier_gates[2 *
                             n] = cell.add_component(f'sublattice_clavier_gates_{2*n}')
            sl_clavier_gates[2*n].component = clavier_gates[2*n]
            sl_clavier_gates[2*n].center = (0, 0)
            sl_clavier_gates[2*n].build()

            clavier_gates[2*n +
                          1] = self.qda_elements.add_clavier_gate(f'clavier_gates_{2*n+1}')
            clavier_gates[2*n+1].layer = clav_layers_order[2*n]
            clavier_gates[2*n+1].width = self.clav_width
            clavier_gates[2*n+1].length = self.clav_length
            clavier_gates[2*n+1].gate_width = self.clav_gate_width
            clavier_gates[2*n+1].gate_length = (self.clav_gate_length[n % 2] +
                                                n*self.clav_gate_width)
            clavier_gates[2*n+1].n_clav_rep = self.n_clav_rep
            clavier_gates[2*n+1].spacing = self.spacing
            clavier_gates[2*n+1].rotation = 180
            clavier_gates[2*n+1].x = ((n+1/2)*self.clav_gate_spacing /
                                      self._n_clav_gates)
            clavier_gates[2*n+1].y = -(self.y+self.clav_gate_length[n % 2]/2 +
                                       self.clav_gap[n % 2]/2 +
                                       n*self.clav_gate_width)
            clavier_gates[2*n+1].fillet = self.fillet
            clavier_gates[2*n+1].fillet_tolerance = self.fillet_tolerance
            clavier_gates[2*n+1].build()

            sl_clavier_gates[2*n +
                             1] = cell.add_component(f'sublattice_clavier_gates_{2*n+1}')
            sl_clavier_gates[2*n+1].component = clavier_gates[2*n+1]
            sl_clavier_gates[2*n+1].center = (0, 0)
            sl_clavier_gates[2*n+1].build()

        screen = self.qda_elements.add_screening_gate('screen_clav')

        screen.layer = self.screen_layer
        screen.vertices = [(-self.screen_length/2, self.screen_width/2),
                           (self.screen_length/2, self.screen_width/2),
                           (self.screen_length/2, -self.screen_width/2),
                           (-self.screen_length/2, -self.screen_width/2)]
        screen.x = 0
        screen.y = -((self.clav_gate_length[0] -
                      self.clav_gap[0] +
                      self.screen_width)/2 +
                     self.screen_gap)
        screen.build()

        sl_clav_screen = cell.add_component(f'sublattice_clav_screen_up')
        sl_clav_screen.component = screen
        sl_clav_screen.center = (0, 0)
        sl_clav_screen.rows = 2
        sl_clav_screen.spacing_y = 2*abs(screen.y)
        sl_clav_screen.build()

        # cell.add(gdstk.Reference(sl_clav_screen.cell))
        self.cell = cell.cell

        return cell


class Sublattice:
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
        self._get_lim(axis=0)
        self._get_lim(axis=1)
        # only update elements attribute if it exists
        if not isinstance(self.component, FanoutGenerator):
            self.elements = update_positions(self.component.elements,
                                             self.rows, self.columns,
                                             self.spacing[0], self.spacing[1],
                                             self.center)
        else:
            pass


class Gate:
    def __init__(self, name):
        self.name = name
        self.points = None
        self.position = None
        self.layer = None


class FanOutLine:
    def __init__(self, name):
        """
        Initialize a FanOutLine object.

        Args:
            name (str): Name of the fanout line.
        """
        self.name = name
        self.element_name = None
        self.n_element = None
        self.layer_fine = None
        self.fo_fine_start = None
        self.fo_fine_end = None
        self.fo_fine_coarse_overlap = None
        self.layer_coarse = None
        self.polygons_fine = None
        self.polygons_coarse = None
        self.fo_direction = None
        self.n_fanout = None
        self.fog = None
        self.cell = gdstk.Cell(self.name)
        self.fillet = 0
        self.fillet_tolerance = 1e-3
        self.elements = {name: {'vertices': [],
                                'positions': [],
                                'layer': self.layer_fine}}

    def assign_fanout(self):
        self.polygons_coarse = self.fog.fo_polygons_coarse[self.fo_direction][self.n_fanout]
        self.fo_fine_start = self.fog.qda_elements[self.element_name]['positions'][self.n_fanout]
        fo_fine_end_poly = self.fog.fo_polygons_coarse[self.fo_direction][self.n_fanout][3:5]
        self.fo_fine_coarse_overlap = self.fog.fo_polygons_coarse[
            self.fo_direction][self.n_fanout][2:6]
        self.fo_fine_end = [int(sum(col) / len(col))
                            for col in zip(*fo_fine_end_poly)]
        self.polygons_fine = self.fo_fine_start

    def build_coarse_fo(self):
        """
        Build the coarse fanout line element.
        """
        fo_line = gdstk.Polygon(self.polygons_coarse, layer=self.layer_coarse)
        fo_line.fillet(self.fillet, tolerance=self.fillet_tolerance)
        fo_line.fillet(0.02, tolerance=1e-4)

        self.elements[self.name]['vertices'] = fo_line.points
        self.elements[self.name]['positions'] = [0, 0]
        self.elements[self.name]['layer'] = self.layer_coarse

        self.cell.add(fo_line)


# %% Fanout Generator

class FanoutGenerator():
    def __init__(self, name, qda):
        # self.layer_dict = None
        # self.fo_width = 0.03
        # self.fo_spacing_north = 0.07
        # self.fo_spacing_east = 0.07
        # self.fo_spacing_south = 0.07
        # self.fo_spacing_west = 0.07
        # self.fo_xoffset_north = 0
        # self.fo_xoffset_south = 0
        # self.fo_yoffset_east = 0
        # self.fo_yoffset_west = 0
        # self.fo2device_buffer = 0.2
        # self.layers_fo_north_dict = {0: list(np.arange(16, dtype=int))}
        # self.layers_fo_east_dict = {0: list(np.arange(16, dtype=int))}
        # self.layers_fo_south_dict = {0: list(np.arange(16, dtype=int))}
        # self.layers_fo_west_dict = {0: list(np.arange(16, dtype=int))}
        # self.chip_width = 4000
        # self.margin = 300
        # self.bondpad_spacing = 100
        # self.bondpad_length = 400
        # self.fan_width = 25
        # self.int_spacing = 1
        # self.int_width = 2
        # self.use_width = self.chip_width - 2 * self.margin - \
        #     2 * self.bondpad_length - self.bondpad_spacing/2
        # self.use_height = self.chip_width - 2 * self.margin
        self.name = name
        self.qda_elements = qda.elements
        self.elements = {}
        self.cell = gdstk.Cell(name)
        self.components = {}
        self.fo_lines = {}
        self._all_directions = ['top', 'bottom', 'right', 'left']
        self.fo_polygons_coarse = None
        self.rect_dims = [(30, 30), (1500, 1500), (2600, 2600)]
        self.fo_widths = [1, 6, 25]
        self.fanout_counts = {'top': 14, 'bottom': 14, 'left': 13, 'right': 13}
        self.spacings = [2, 40, 160]
        self.fo_fine_coarse_overlap = 3
        self.bondpad_position = {'top': 3000, 'bottom': 3000,
                                 'left': 3000, 'right': 3000}
        self.bondpad_size = {'top': (110, 400), 'bottom': (110, 400),
                             'left': (400, 110), 'right': (400, 110)}
        self._n = 0

    def create_fo_polygons_coarse(self):
        polygons = {}
        fanout_positions = compute_fanout_positions(self.rect_dims,
                                                    self.fanout_counts,
                                                    self.spacings)
        fo_lines = get_fo_lines(fanout_positions,
                                self.fanout_counts, self.rect_dims)
        for direction in self._all_directions:
            polygons[direction] = [generate_polygon_for_fanout(direction,
                                                               fo_lines[direction][n_fo],
                                                               self.bondpad_position,
                                                               self.bondpad_size,
                                                               self.fo_widths,
                                                               self.fo_fine_coarse_overlap)
                                   for n_fo in range(self.fanout_counts[direction])]

        self.fo_polygons_coarse = polygons

    # def create_fine_fo_polygons(self):

    # def add_fanout_line(self, name):
    #     """
    #     Add a fanout element to the collection.

    #     Args:
    #         name: Name of the fanout element.

    #     Returns:
    #         The created fanout element.

    #     """
    #     fanout_line = FanOutLine(name)
    #     self.components[name] = fanout_line
    #     return fanout_line

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
            self.cell.add(gdstk.Reference(cell.component.cell))
            elements = merge_dicts(elements, cell.component.elements)
        self.elements = elements
        self.cell.flatten()
