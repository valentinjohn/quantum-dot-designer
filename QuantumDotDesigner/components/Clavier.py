# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:46:38 2023

@author: vjohn
"""
# %% imports

from QuantumDotDesigner.base import UnitCell
from QuantumDotDesigner.BaseCollection import BaseCollection

from QuantumDotDesigner.elements.ClavierGate import ClavierGate
from QuantumDotDesigner.elements.ScreeningGate import ScreeningGate

import numpy as np

# %% definition


class Clavier(UnitCell):
    def __init__(self, name, collection: BaseCollection):
        # if not isinstance(qda_elements, QuantumDotArrayElements):
        #     raise TypeError(
        #         f"Expected qda_elements to be of type {QuantumDotArrayElements}, but got {type(qda_elements)} instead.")
        super().__init__(name)
        self.collection = collection
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

        self.collection.add_component(self)

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

        self.clavier_gates[name] = ClavierGate(name, self.collection)
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

        self.clavier_gates[name_odd] = ClavierGate(name_odd, self.collection)
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
        self.screen = ScreeningGate(f'{self.name}_screen', self.collection)
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
