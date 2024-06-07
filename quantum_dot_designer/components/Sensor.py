# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:33:18 2023

@author: vjohn
"""

# %% imports

from ..base import Component
from ..BaseCollection import BaseCollection

from ..elements.Plunger import Plunger
from ..elements.Barrier import Barrier
from ..elements.Ohmic import Ohmic



import numpy as np
import copy

# %% definition


class Sensor(Component):
    """
    Represents a sensor component in a quantum dot design, containing elements like plungers, barriers, and ohmics.

    The `Sensor` class is an extension of the `Component` class and represents a more complex entity in the quantum dot system. It manages the geometric and spatial configuration of its constituent elements and handles their collective behavior as a single component.

    Attributes:
        gap_ohmic_pl (int): Gap between the ohmic and plunger elements.
        gap_sep (int): Gap between the separation gate and plunger elements.
        source_pos (str): Position identifier for the source element (e.g., 'top', 'bottom').
        source_pos_angle (float): Angular position of the source element, in radians.
        drain_pos (str): Position identifier for the drain element.
        drain_pos_angle (float): Angular position of the drain element, in radians.
        sep_pos (str): Position identifier for the separation gate.
        sep_pos_angle (float): Angular position of the separation gate, in radians.
        barrier_orientation (dict): Orientation specifications for barriers.
        offset (dict): Offset values for various elements within the sensor.
        element_positions (dict): Stores the positions of constituent elements.
        _bar_angle_dict (dict): Dictionary mapping position identifiers to angles.
        components_position (dict): Positions of the various components within the sensor.
        fillet (tuple): Parameters for the fillet operation in geometry construction.
        __feature_gap (float): Internal attribute to manage feature gaps.

    Methods:
        _init_elements(name, collection): Initializes the elements that make up the sensor.
        _calculate_positions(): Calculates the positions of the elements based on sensor geometry.
        _set_barrier_properties(): Sets properties related to barriers within the sensor.
        _build_and_add_elements(): Builds and adds elements to the sensor component.
        _build_elements(): Handles the comprehensive building process for the sensor elements.
        copy(copy_name, collection): Creates a copy of the sensor with a new name and potentially new attributes.
        build(): Overridden method that builds the sensor component, assembling its subcomponents and finalizing its geometry.
    """

    def __init__(self, name: str, collection: BaseCollection):
        if not isinstance(collection, BaseCollection):
            raise TypeError(
                f"Expected collection to be of type {BaseCollection}, but got {type(collection)} instead.")

        super().__init__(name)
        self._init_elements(name, collection)

        self.gap_ohmic_pl = 40
        self.gap_sep = 40
        self.source_pos = None
        self.source_pos_angle = None
        self.drain_pos = None
        self.drain_pos_angle = None
        self.sep_pos = None
        self.sep_pos_angle = None
        self.barrier_orientation = {'source': 'counter-clockwise',
                                    'drain': 'clockwise',
                                    'sep': 'clockwise',
                                    }
        self.offset = {'source': (0, 0),
                       'drain': (0, 0),
                       'barrier_source': (0, 0),
                       'barrier_drain': (0, 0),
                       }
        self.element_positions = {'source': None,
                                  'drain': None,
                                  'barrier_source': None,
                                  'barrier_drain': None,
                                  'barrier_sep': None
                                  }
        self._bar_angle_dict = {'top': 0, 'right': np.pi/2,
                                'bottom': np.pi, 'left': -np.pi/2,
                                'top-right': np.pi/4,
                                'bottom-right': 3/4*np.pi,
                                'bottom-left': -3/4*np.pi,
                                'top-left': -1/4*np.pi}
        self.components_position = {}
        self.fillet = (0, 1e-3)
        self.__feature_gap = None
        collection.add_component(self)

    def _init_elements(self, name, collection):
        self.plunger = Plunger(f'{name}_plunger', collection)
        self.barrier_source = Barrier(f'{name}_barrier_source', collection)
        self.barrier_drain = Barrier(f'{name}_barrier_drain', collection)
        self.source = Ohmic(f'{name}_source', collection)
        self.drain = Ohmic(f'{name}_drain', collection)
        self.barrier_sep = Barrier(f'{name}_barrier_seperation', collection)

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

        self.element_positions['source'] = ((i*(plunger._asymx*plunger.diameter/2 +
                                                self.gap_ohmic_pl) +
                                             self.offset['source'][0],
                                             j*(plunger._asymy*plunger.diameter/2 +
                                                self.gap_ohmic_pl) +
                                             self.offset['source'][1]))
        self.element_positions['drain'] = ((m*(plunger._asymx*plunger.diameter/2 +
                                               self.gap_ohmic_pl) +
                                            self.offset['drain'][0],
                                            n*(plunger._asymy*plunger.diameter/2 +
                                               self.gap_ohmic_pl) +
                                            self.offset['drain'][1]))

        self.element_positions['barrier_source'] = (i*(plunger._asymx*plunger.diameter/2 +
                                                       bar_source.width/2-self.__feature_gap) +
                                                    self.offset['barrier_source'][0],
                                                    j*(plunger._asymy*plunger.diameter/2 +
                                                       bar_source.width/2-self.__feature_gap) +
                                                    self.offset['barrier_source'][1])
        self.element_positions['barrier_drain'] = (m*(plunger._asymx*plunger.diameter/2 +
                                                      bar_drain.width/2-self.__feature_gap) +
                                                   self.offset['barrier_drain'][0],
                                                   n*(plunger._asymy*plunger.diameter/2 +
                                                      bar_drain.width/2 - self.__feature_gap) +
                                                   self.offset['barrier_drain'][1])

        self.element_positions['barrier_sep'] = (u*(plunger._asymx*plunger.diameter/2+self.gap_sep/2),
                                                 v*(plunger._asymy*plunger.diameter/2+self.gap_sep/2))

    def _set_barrier_properties(self):
        source = self.source
        drain = self.drain
        bar_source = self.barrier_source
        bar_drain = self.barrier_drain
        bar_sep = self.barrier_sep

        if self.barrier_orientation['drain'] == 'clockwise':
            drain.sensor_pos = 'top'
            drain.ohmic_pos = 'right'
        else:
            drain.sensor_pos = 'bottom'
            drain.ohmic_pos = 'right'

        if self.barrier_orientation['source'] == 'clockwise':
            source.sensor_pos = 'top'
            source.ohmic_pos = 'right'
        else:
            source.sensor_pos = 'bottom'
            source.ohmic_pos = 'right'

        bar_drain_angle_offset = 0
        multiplier_bar_drain = 1
        if self.barrier_orientation['drain'] == 'clockwise':
            bar_drain_angle_offset = np.pi
            multiplier_bar_drain = -1

        bar_source_angle_offset = 0
        multiplier_bar_source = 1
        if self.barrier_orientation['source'] == 'clockwise':
            bar_source_angle_offset = np.pi
            multiplier_bar_source = -1

        bar_sep_angle_offset = 0
        if self.barrier_orientation['sep'] == 'clockwise':
            bar_sep_angle_offset = np.pi

        bar_source.rotate = -self.source_pos_angle+bar_source_angle_offset
        source.rotate = (-self.source_pos_angle+bar_source_angle_offset +
                         multiplier_bar_source*np.pi/2)

        bar_drain.rotate = -self.drain_pos_angle+bar_drain_angle_offset
        drain.rotate = (-self.drain_pos_angle+bar_drain_angle_offset +
                        multiplier_bar_drain*np.pi/2)

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
        sl_bs.center = self.element_positions['barrier_source']
        sl_bs.build()

        sl_bd = self.add_component(self.barrier_drain)
        sl_bd.center = self.element_positions['barrier_drain']
        sl_bd.build()

        sl_bsep = self.add_component(self.barrier_sep)
        sl_bsep.center = self.element_positions['barrier_sep']
        sl_bsep.build()

        sl_source = self.add_component(self.source)
        sl_source.center = self.element_positions['source']
        sl_source.build()

        sl_drain = self.add_component(self.drain)
        sl_drain.center = self.element_positions['drain']
        sl_drain.build()

    def _build_elements(self):
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

    def copy(self, copy_name, collection):
        # if not component.built:
        #     component.build()
        attributes = copy.copy(vars(self))
        attributes.pop('name')
        attributes.pop('cell')
        attributes.pop('elements')
        attributes.pop('components')

        new_element = type(self)(copy_name, collection)

        new_element.__dict__.update(attributes)

        new_element.cell = new_element.cell.copy(copy_name)
        new_element.plunger = self.plunger.copy(f'{copy_name}_plunger',
                                                collection)
        new_element.barrier_source = self.barrier_source.copy(f'{copy_name}_barrier_source',
                                                              collection)
        new_element.barrier_drain = self.barrier_drain.copy(f'{copy_name}_barrier_drain',
                                                            collection)
        new_element.source = self.source.copy(
            f'{copy_name}_source', collection)
        new_element.drain = self.drain.copy(f'{copy_name}_drain', collection)
        new_element.barrier_sep = self.barrier_sep.copy(f'{copy_name}_barrier_seperation',
                                                        collection)

        return new_element

    def build(self):
        self._build_elements()
        super().build()
