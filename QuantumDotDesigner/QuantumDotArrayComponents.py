# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 18:02:07 2023

@author: vjohn
"""

import copy

from QuantumDotDesigner import QuantumDotArrayElements
from QuantumDotDesigner.components.Sensor import Sensor
from QuantumDotDesigner.components.Clavier import Clavier
from QuantumDotDesigner.components.FanOutLine import FanOutLine


class QuantumDotArrayComponents:
    def __init__(self, qda_elements):
        if not isinstance(qda_elements, QuantumDotArrayElements):
            raise TypeError(
                f"Expected qda to be of type {QuantumDotArrayElements}, but got {type(qda_elements)} instead.")
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

    def add_fo_line(self, element_name, element_number=0):
        fanout_line = FanOutLine(
            element_name, element_number, self.qda_elements)
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
            new_element.source = qda_elements.add_copy(component.source,
                                                       f'{copy_name}_source')
            new_element.drain = qda_elements.add_copy(component.drain,
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
