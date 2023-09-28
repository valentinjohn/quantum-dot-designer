# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:07:21 2023

@author: vjohn
"""

# %% imports

from abc import abstractmethod
from QuantumDotDesigner.base.ElementBase import ElementBase
from QuantumDotDesigner.BaseCollection import BaseCollection

import copy

# %% definition


class Element(ElementBase):
    def __init__(self, name: str, collection: BaseCollection):
        """
        Initialize an Element object.

        Args:
            name (str): Name of the element.
            layer (int): Description of layer. Default is None.
            rotate (float): Rotation of the element. Default is 0.0.
            fillet (float): Fillet value. Indicates how much the polygon is rounded.
        """

        super().__init__(name)
        collection.add_element(self)

    def copy(self, copy_name, collection: BaseCollection):
        # if not component.built:
        #     component.build()
        attributes = copy.copy(vars(self))
        attributes.pop('name')
        attributes.pop('cell')
        attributes.pop('elements')
        # attributes.pop('components')

        new_element = type(self)(copy_name, collection)
        new_element.__dict__.update(attributes)
        collection.add_element(new_element)

        return new_element
