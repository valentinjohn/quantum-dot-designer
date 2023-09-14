# -*- coding: utf-8 -*-
"""
Created on Tue May  9 07:55:18 2023

"""

# %% imports

import gdstk
import numpy as np
import copy

from QuantumDotDesigner.helpers.helpers import *

from QuantumDotDesigner.base.UnitCell import UnitCell

from QuantumDotDesigner.elements.Plunger import Plunger
from QuantumDotDesigner.elements.Barrier import Barrier
from QuantumDotDesigner.elements.ScreeningGate import ScreeningGate
from QuantumDotDesigner.elements.Ohmic import Ohmic
from QuantumDotDesigner.elements.ClavierGate import ClavierGate
from QuantumDotDesigner.elements.FanOutLineFine import FanOutLineFine
from QuantumDotDesigner.elements.FanOutLineCoarse import FanOutLineCoarse

# %% Quantum Dot Array Elements


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

    # def add_sensor(self, name):
    #     """
    #     Add a sensor to the collection.

    #     Args:
    #         name: Name of the sensor.

    #     Returns:
    #         The created sensor.

    #     """
    #     sensor = Sensor(name, self)
    #     self.components[name] = sensor
    #     return sensor

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

    # def add_clavier(self, name):
    #     """
    #     Add a clavier to the collection.

    #     Args:
    #         name: Name of the clavier.

    #     Returns:
    #         The created clavier.

    #     """
    #     clavier = Clavier(name, self)
    #     self.components[name] = clavier
    #     return clavier

    def add_fo_line_fine(self, name):
        fanout_line = FanOutLineFine(name)
        self.components[name] = fanout_line
        return fanout_line

    def add_fo_line_coarse(self, name):
        fanout_line = FanOutLineCoarse(name)
        self.components[name] = fanout_line
        return fanout_line

    # def add_fo_line(self, element_name, element_number):
    #     fanout_line = FanOutLine(element_name, element_number, self)
    #     self.components[f'fo_{element_name}_{element_number}'] = fanout_line
    #     return fanout_line

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


class QuantumDotArrayComponents:

    def __init__(self, qda_elements):
        self.components = {}
        self.qda_elements = qda_elements

    def add_sensor(self, name):
        """
        Add a sensor to the collection.

        Args:
            name: Name of the sensor.

        Returns:
            The created sensor.

        """
        sensor = Sensor(name, self.qda_elements)
        self.components[name] = sensor
        return sensor

    def add_clavier(self, name):
        """
        Add a clavier to the collection.

        Args:
            name: Name of the clavier.

        Returns:
            The created clavier.

        """
        clavier = Clavier(name, self.qda_elements)
        self.components[name] = clavier
        return clavier

    def add_fo_line(self, element_name, element_number):
        fanout_line = FanOutLine(
            element_name, element_number, self.qda_elements)
        self.components[f'fo_{element_name}_{element_number}'] = fanout_line
        return fanout_line

    def add_copy(self, component, copy_name):
        if not component.built:
            component.build()
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


# %% Components


class Sensor(UnitCell):
    def __init__(self, name: str, qda_elements):
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
        self.fo_points = None
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
        self.elements[self.fo_line_fine.name] = {'vertices': [],
                                                 'positions': [],
                                                 'layer': self.fo_line_fine.layer}

        self.fo_line_fine.fo_direction = self.fo_direction
        self.fo_line_fine.layer = self.qda_elements.components[self.element_name].layer

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
