# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:07:21 2023

@author: vjohn
"""

# %% imports

from .ElementBase import ElementBase
from ..BaseCollection import BaseCollection


import copy

# %% definition


class Element(ElementBase):
    """
    Represents a specific type of element in the quantum_dot_designer system, extending the base functionalities
    provided by the ElementBase class. Each Element object is registered within a collection upon initialization.

    The Element class integrates with the BaseCollection, allowing each element to be part of a broader collection 
    for organized storage and management. It also overrides the copy functionality to handle deep copying of 
    the element and its attributes, ensuring the new element is also registered within the desired collection.

    Attributes:
        Inherits all attributes from ElementBase.

    Methods:
        __init__(self, name: str, collection: BaseCollection): Initializes the Element and registers it within the provided collection.
        copy(self, copy_name: str, collection: BaseCollection): Creates a deep copy of the element with a new name and registers it within the provided collection.
    """

    def __init__(self, name: str, collection: BaseCollection):
        """
        Initialize an Element object and register it to a collection.

        The constructor takes in the name of the element and a collection to which the element will belong. 
        It calls the superclass constructor to set up the element itself, then adds the element to the 
        specified collection.

        Args:
            name (str): Unique identifier for the element.
            collection (BaseCollection): The collection to which the element will be added.

        """
        super().__init__(name)
        collection.add_element(self)

    def copy(self, copy_name, collection: BaseCollection):
        attributes = copy.copy(vars(self))
        attributes.pop('name')
        attributes.pop('cell')
        attributes.pop('elements')

        new_element = type(self)(copy_name, collection)
        new_element.__dict__.update(attributes)
        collection.add_element(new_element)

        return new_element
