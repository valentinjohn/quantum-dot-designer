# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 08:38:11 2023

@author: vjohn
"""

from .base.ElementBase import ElementBase
from .base.UnitCell import UnitCell
from .base.Layer import Layer


class BaseCollection:
    """
A storage class that maintains collections of elements, components, and layers for easy access and management.

This class holds dictionaries of different entities and provides methods to add elements, components, 
and layers to their respective dictionaries. Each entity is stored in a dictionary with its name as 
the key for quick retrieval.

Attributes:
    elements (dict): A dictionary of elements, with element names as keys and corresponding element objects as values.
    components (dict): A dictionary of components, with component names as keys and corresponding component objects as values.
    fo_elements (dict): A dictionary of fanout elements, with their names as keys and corresponding objects as values.
    fo_components (dict): A dictionary of fanout components, with their names as keys and corresponding objects as values.
    layers (dict): A dictionary of layers, with layer names as keys and corresponding layer objects as values.

Methods:
    add_element(element: ElementBase): Adds an element to the 'elements' dictionary.
    add_component(component: UnitCell): Adds a component to the 'components' dictionary.
    add_fo_element(element: ElementBase): Adds an 'fo_element' to the 'fo_elements' dictionary.
    add_fo_component(component: UnitCell): Adds an 'fo_component' to the 'fo_components' dictionary.
    add_layer(layer: Layer): Adds a layer to the 'layers' dictionary.
"""

    def __init__(self):
        self.elements = {}
        self.components = {}
        self.fo_elements = {}
        self.fo_components = {}
        self.layers = {}

    def add_element(self, element: ElementBase):
        self.elements[element.name] = element

    def add_component(self, component: UnitCell):
        self.components[component.name] = component

    def add_fo_element(self, element: ElementBase):
        self.fo_elements[element.name] = element

    def add_fo_component(self, component: UnitCell):
        self.fo_components[component.name] = component

    def add_layer(self, layer: Layer):
        self.layers[layer.name] = layer
