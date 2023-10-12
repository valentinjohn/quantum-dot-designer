# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:07:21 2023

@author: vjohn
"""

# %% imports

from QuantumDotDesigner import BaseCollection
import copy

# %% definition


class Layer:
    def __init__(self, name: str, collection: BaseCollection,
                 fine=0, coarse=1, via_etch=None, via_fine=None, via_coarse=None):
        self.name = name
        self.fine = fine
        self.coarse = coarse
        self.via_etch = via_etch
        self.via_fine = via_fine
        self.via_coarse = via_coarse
        collection.add_layer(self)

    def copy(self, copy_name, collection: BaseCollection):
        attributes = copy.copy(vars(self))
        attributes.pop('name')

        new_element = type(self)(copy_name, collection)
        new_element.__dict__.update(attributes)
        collection.add_element(new_element)

        return new_element
