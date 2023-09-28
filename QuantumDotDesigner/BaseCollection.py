# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 08:38:11 2023

@author: vjohn
"""

from QuantumDotDesigner.base.ElementBase import ElementBase
from QuantumDotDesigner.base.UnitCell import UnitCell


class BaseCollection:
    def __init__(self):
        self.elements = {}
        self.components = {}

    def add_element(self, element: ElementBase):
        self.elements[element.name] = element

    def add_component(self, component: UnitCell):
        self.components[component.name] = component
