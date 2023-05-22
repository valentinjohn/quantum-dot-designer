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

#%% imports

import gdstk
import numpy as np
import copy
from abc import ABC, abstractmethod


#%% definitions

def rot_mat(theta):
    """
    Compute the rotation matrix for a given angle.

    Args:
        theta: Angle of rotation in radians.

    Returns:
        Rotation matrix as a numpy array.

    """
    return np.array([[np.cos(theta),-np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])

def gen_poly(n, sp = None):
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
        self.elements = {}
        self.sublattices = []
        self.unit_cells = {}

    def add_plunger(self, name):
        """
        Add a plunger element to the collection.

        Args:
            name: Name of the plunger element.

        Returns:
            The created plunger element.

        """
        plunger = Plunger(name)
        self.elements[name] = plunger
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
        self.elements[name] = barrier
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
        self.elements[name] = screening_gate
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
        self.elements[name] = ohmic
        return ohmic

    def add_sensor(self, name):
        """
        Add a sensor to the collection.

        Args:
            name: Name of the sensor.

        Returns:
            The created sensor.

        """
        sensor = Sensor(name)
        self.elements[name] = sensor
        return sensor


class UnitCell:
    def __init__(self, qda_instance, name='unit_cell'):
        """
        Initialize a UnitCell object.

        Args:
            qda_instance: The QuantumDotArray instance that this UnitCell belongs to.
            name (str): Name of the unit cell (default: 'unit_cell').
        """
        self.parent = qda_instance
        self.sublattices = {}
        self.elements = {}
        self.cells = []
        self.unit_cell = gdstk.Cell(name)

    def add_sublattice(self, name):
        """
        Add a sublattice to the unit cell.

        Args:
            name (str): Name of the sublattice.

        Returns:
            Sublattice: The created Sublattice object.
        """
        sublattice = Sublattice(name)
        self.sublattices[name] = sublattice
        self.cells.append(sublattice)
        return sublattice

    def build(self):
        """
        Build the unit cell by adding sublattices to the unit cell.

        This method adds the sublattices to the unit cell's cell object.
        """
        for cell in self.cells:
            self.unit_cell.add(gdstk.Reference(cell.cell))
        self.unit_cell.flatten()


class QuantumDotArray:
    def __init__(self, parent_instance, unitcell_instance=None):
        """
        Initialize a QuantumDotArray object.

        Args:
            parent_instance: The parent instance.
            unitcell_instance (UnitCell): An optional UnitCell instance to use for the array.
        """
        self.spacing_qd = 200
        self.spacing_qd_diag = 2**0.5 * self.spacing_qd
        self.parent_instance = parent_instance
        self.unitcell_instance = unitcell_instance
        self.sublattices = {}
        self.cells = []
        self.elements = self.parent_instance.elements
        self.main_cell = gdstk.Cell('MAIN')
        self.add_unitcells()

    def add_unitcells(self, name='unit_cell'):
        """
        Add unit cells to the QuantumDotArray.

        Args:
            name (str): Name of the unit cells (default: 'unit_cell').
        """
        for cell in self.unitcell_instance:
            self.cells.append(cell)

    def add_sublattice(self, name):
        """
        Add a sublattice to the QuantumDotArray.

        Args:
            name (str): Name of the sublattice.

        Returns:
            Sublattice: The created Sublattice object.
        """
        sublattice = Sublattice(name)
        self.sublattices[name] = sublattice
        self.cells.append(sublattice)
        return sublattice

    def build(self):
        """
        Build the QuantumDotArray.

        This method adds the sublattices and unit cells to the main cell object.
        """
        for cell in self.cells:
            if isinstance(cell, Sublattice):
                self.main_cell.add(gdstk.Reference(cell.cell))
            elif isinstance(cell, UnitCell):
                for c in cell.cells:
                    self.main_cell.add(gdstk.Reference(c.cell))
        self.main_cell.flatten()

    def save_as_gds(self, filename):
        """
        Save the QuantumDotArray as a GDS file.

        Args:
            filename (str): Name of the output GDS file.
        """
        lib = gdstk.Library()
        lib.add(self.main_cell)
        lib.write_gds(filename)


class Element(ABC):
    def __init__(self, name):
        """
        Initialize an Element object.

        Args:
            name (str): Name of the element.
        """
        self.name = name
        self.layer = None
        self.polygons = []
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

    def copy(self, copy_name):
        """
        Create a copy of the element.

        Args:
            copy_name (str): Name of the copied element.

        Returns:
            Element: The copied Element object.
        """
        attributes = copy.copy(vars(self))
        attributes.pop('name')
        attributes.pop('cell')
        new_element = type(self)(copy_name)
        new_element.__dict__.update(attributes)
        return new_element


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

    def update_asym(self, asym):
        """
        Update the asymmetry of the plunger.

        Args:
            asym (float): Asymmetry value.
        """
        self._asymx = asym
        self._asymy = 1 / asym

    def get_asym(self):
        """
        Get the asymmetry of the plunger.

        Returns:
            tuple: Tuple containing the X and Y asymmetry values.
        """
        self.update_asym(self.asym)
        return (self._asymx,  self._asymy)

    def build(self):
        """
        Build the plunger element.
        """
        self.update_asym(self.asym)
        pl_points = gen_poly(8)
        pl = gdstk.Polygon(pl_points, layer=self.layer)
        pl.scale(0.5 / np.cos(np.pi / 8) * self.diameter)
        pl.scale(sx=self._asymx, sy=self._asymy)
        pl.translate(self.x, self.y)
        pl.fillet(self.fillet, tolerance=self.fillet_tolerance)
        pl.fillet(0.02, tolerance=1e-4)
        cell = gdstk.Cell(self.name)
        cell.add(pl)
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
        self.shape = None

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
        self.cell = cell


class ScreeningGate(Element):
    def __init__(self, name):
        super().__init__(name)
        self.shape = None

    def build(self):
        print('Build method for Screening gate not implemeneted yet.')

class Ohmic(Element):
    def __init__(self, name):
        super().__init__(name)
        self.width = 0.0
        self.height = 0.0
        self.shape = None

    def build(self):
        print('Build method for Ohmic not implemeneted yet.')

class Sensor:
    def __init__(self, name):
        """
        Initialize a Sensor object.

        Args:
            name (str): Name of the sensor.
        """
        self.name = name
        self.cell = None
        self.x = 0
        self.y = 0
        self.plunger = Plunger(f'{name} plunger')
        self.barrier_source = Barrier(f'{name} barrier source')
        self.barrier_drain = Barrier(f'{name} barrier drain')
        self.source = Ohmic(f'{name} source')
        self.drain = Ohmic(f'{name} barrier drain')
        self.barrier_sep = Barrier(f'{name} barrier seperation')
        self.gap_ohmic_pl = 40
        self.gap_sep = 40
        self.source_pos = 'left'
        self.drain_pos = 'right'
        self.sep_pos = 'bottom'
        self.source_position_offset = (0,0)
        self.drain_position_offset = (0,0)
        self.bar_sou_position_offset = (0,0)
        self.bar_dra_position_offset = (0,0)
        self.bar_sharp_source = (0,0)
        self.bar_sharp_drain = (0,0)
        self.fillet = (0, 1e-3)
        self.__feature_gap = None
        self.sd_position = None
        self.bar_position = None

    def copy(self, copy_name):
        """
        Create a copy of the sensor.

        Args:
            copy_name (str): Name of the copied sensor.

        Returns:
            Sensor: The copied Sensor object.
        """
        attributes = copy.copy(vars(self))
        attributes.pop('name')
        attributes.pop('cell')
        new_element = type(self)(copy_name)
        new_element.__dict__.update(attributes)
        return new_element

    def build(self):
        """
        Build the sensor element.

        Returns:
            gdstk.Cell: The built sensor cell.
        """
        cell = gdstk.Cell(self.name)
        plunger = self.plunger
        bar_source = self.barrier_source
        bar_drain = self.barrier_drain
        source = self.source
        drain = self.drain
        bar_sep = self.barrier_sep
        self.__feature_gap = self.barrier_source.width - self.gap_ohmic_pl

        orientation_dict = {'top':(0,1), 'right':(1,0),
                            'bottom':(0,-1), 'left':(-1,0),
                            'top-right':(2**0.5/2,2**0.5/2),
                            'bottom-right':(2**0.5/2,-2**0.5/2),
                            'bottom-left':(-2**0.5/2,-2**0.5/2),
                            'top-left':(-2**0.5/2,2**0.5/2)}

        bar_angle_dict = {'top':np.pi, 'right':-np.pi/2,
                          'bottom':np.pi, 'left':+np.pi/2,
                          'top-right':-np.pi/4,
                          'bottom-right':np.pi/4,
                          'bottom-left':3*np.pi/4,
                          'top-left':-3*np.pi/4}

        (i,j) = orientation_dict[self.source_pos]
        (m,n) = orientation_dict[self.drain_pos]
        (u,v) = orientation_dict[self.sep_pos]

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

        self.plunger.build()
        self.barrier_source.build()
        self.barrier_drain.build()
        # self.source.build()
        # self.drain.build()
        self.barrier_sep.build()

        cell.add(gdstk.Reference(plunger.cell))
        cell.add(gdstk.Reference(bar_source.cell))
        cell.add(gdstk.Reference(bar_drain.cell))
        # cell.add(gdstk.Reference(source))
        # cell.add(gdstk.Reference(drain))
        cell.add(gdstk.Reference(bar_sep.cell))
        self.cell = cell

        return cell


class Sublattice:
    def __init__(self, name):
        """
        Initialize a Sublattice object.

        Args:
            name (str): Name of the sublattice.
        """
        self.name = name
        self.element = None
        self.n_rows = 1
        self.n_columns = 1
        self.spacing_x = 100
        self.spacing_y = 100
        self.center = (0, 0)
        self.positions = list()
        self.points = list()
        self.cell = None
        self._width = (self.n_columns-1) * self.spacing_x
        self._height = (self.n_rows-1) * self.spacing_y
        self.xmax = None
        self.ymax = None

    def _update_width(self):
        """
        Update the width of the sublattice based on the number of columns and spacing.
        """
        self._width = (self.n_columns-1) * self.spacing_x

    def _update_height(self):
        """
        Update the height of the sublattice based on the number of rows and spacing.
        """
        self._height = (self.n_rows-1) * self.spacing_y

    def set_element(self, element):
        """
        Set the element of the sublattice.

        Args:
            element (Element): The element to be placed in the sublattice.
        """
        self.element = element

    def set_rows(self, rows):
        """
        Set the number of rows in the sublattice.

        Args:
            rows (int): Number of rows.
        """
        self.n_rows = rows
        self._update_height()

    def set_columns(self, columns):
        """
        Set the number of columns in the sublattice.

        Args:
            columns (int): Number of columns.
        """
        self.n_columns = columns
        self._update_width()

    def set_center(self, center):
        """
        Set the center position of the sublattice.

        Args:
            center (tuple): Center position as a tuple (x, y).
        """
        self.center = center

    def set_xspacing(self, spacing_x):
        """
        Set the horizontal spacing between elements in the sublattice.

        Args:
            spacing_x (float): Horizontal spacing.
        """
        self.spacing_x = spacing_x
        self._update_width()

    def set_yspacing(self, spacing_y):
        """
        Set the vertical spacing between elements in the sublattice.

        Args:
            spacing_y (float): Vertical spacing.
        """
        self.spacing_y = spacing_y
        self._update_height()

    def get_lim(self, axis=0):
        """
        Get the minimum and maximum limits along the specified axis of the sublattice cell.

        Args:
            axis (int): Axis along which to find the limits (0 for x-axis, 1 for y-axis).

        Returns:
            tuple: Minimum and maximum limits as a tuple (min, max).
        """
        cell_max = 0
        cell_min = 0
        for poly in self.cell.get_polygons():
            poly_max = poly.points[:,axis].max()
            poly_min = poly.points[:,axis].min()
            cell_max = max(cell_max, poly_max)
            cell_min = min(cell_min, poly_min)
        return (cell_min, cell_max)

    def get_positions(self):
        """
        Get the positions of the elements in the sublattice.

        Returns:
            list: List of element positions as tuples (x, y).
        """
        positions = list()
        x0 = self.center[0] - self._width/2
        y0 = self.center[1] - self._height/2
        for row in reversed(range(self.n_rows)):
            for col in range(self.n_columns):
                x = x0 + col*self.spacing_x
                y = y0 + row*self.spacing_y
                positions.append([x, y])
        self.positions = positions

        return positions

    def get_points(self):
        """
        Get the points defining the elements in the sublattice.

        Returns:
            list: List of polygons points defining the elements in the sublattice.
        """
        poly_points = list()
        x0 = self.center[0] - self._width/2
        y0 = self.center[1] - self._height/2
        polygons = self.element.cell.polygons
        for row in reversed(range(self.n_rows)):
            for col in range(self.n_columns):
                for poly in polygons:
                    x = x0 + col*self.spacing_x
                    y = y0 + row*self.spacing_y
                    poly_points.append(poly.points + [x, y])

        self.points = poly_points

        return poly_points

    def build(self):
        """
        Build the sublattice cell by adding the element references to the sublattice cell.
        """
        self.get_positions()
        self.get_points()
        cell = gdstk.Cell(self.name)
        cell.add(gdstk.Reference(self.element.cell,
                                 (self.center[0] - self._width/2,
                                  self.center[1] - self._height/2),
                                 columns=self.n_columns,
                                 rows=self.n_rows,
                                 spacing=(self.spacing_x, self.spacing_y)
                                 ))
        self.cell = cell
        self.xlim = self.get_lim(axis=0)
        self.ylim = self.get_lim(axis=1)


class Gate:
    def __init__(self, name):
        self.name = name
        self.points = None
        self.position = None
        self.layer = None
