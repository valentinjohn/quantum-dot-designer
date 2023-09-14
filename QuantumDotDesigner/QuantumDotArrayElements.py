# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 17:49:40 2023

@author: vjohn
"""

import copy

from QuantumDotDesigner.elements.Plunger import Plunger
from QuantumDotDesigner.elements.Barrier import Barrier
from QuantumDotDesigner.elements.ScreeningGate import ScreeningGate
from QuantumDotDesigner.elements.Ohmic import Ohmic
from QuantumDotDesigner.elements.ClavierGate import ClavierGate
from QuantumDotDesigner.elements.FanOutLineFine import FanOutLineFine
from QuantumDotDesigner.elements.FanOutLineCoarse import FanOutLineCoarse


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

    def add_fo_line_fine(self, name):
        fanout_line = FanOutLineFine(name)
        self.components[name] = fanout_line
        return fanout_line

    def add_fo_line_coarse(self, name):
        fanout_line = FanOutLineCoarse(name)
        self.components[name] = fanout_line
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

        self.components[copy_name] = new_element
        return new_element
