# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 08:38:11 2023

@author: vjohn
"""

from QuantumDotDesigner.base import ElementBase, UnitCell, Layer


class BaseCollection:
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
