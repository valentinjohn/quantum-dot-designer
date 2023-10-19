# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from QuantumDotDesigner.base import Element
from QuantumDotDesigner.helpers.helpers import generate_clavier_gates
from QuantumDotDesigner.BaseCollection import BaseCollection
import gdstk


class ClavierGate(Element):
    """
    Represents a Clavier gate within the QuantumDotDesigner system, a specialized structure used in quantum dot device designs.

    The ClavierGate class extends the standard Element class, introducing specific parameters that define the unique geometry of a Clavier gate. This gate features a series of repeating structures resembling a keyboard, which can be used for conveyor-mode shuttling of electrons and holes.
    This class allows for detailed customization of the gate's properties, including dimensions, spacing, and repetition count. The build process involves generating the gate's geometric representation and incorporating it into the device design.

    Attributes:
        layer (Layer): The layer to which the Clavier gate belongs. It must be a valid Layer instance.
        layer_stage (str): The stage of the layer associated with the Clavier gate, defaults to 'fine'.
        width (int): The width of the individual gate segments.
        length (int): The total length of the Clavier gate structure.
        gate_width (int): The width of the gate structure.
        gate_length (int): The length of the individual gate segments.
        n_clav_rep (int): The number of repeating fingers in the Clavier gate structure.
        spacing (int): The distance between individual fingers in the Clavier gate.
        shift (int): Lateral displacement of the Clavier gate structure.
        x (int): The x-coordinate of the Clavier gate's position.
        y (int): The y-coordinate of the Clavier gate's position.
        fillet (float): The radius of the corners of the Clavier gate's structure.
        fillet_tolerance (float): The precision used in constructing the filleted corners.
        rotation (float): The angle of rotation applied to the Clavier gate structure.

    Methods:
        build(): Constructs the Clavier gate's geometric representation and integrates it into the device design.
    """

    def __init__(self, name, collection: BaseCollection):
        """
        Initialize a Clavier gate object.

        Args:
            name (str): Name of the clavier gate
        """
        super().__init__(name, collection)
        self.layer = None
        self.layer_stage = 'fine'
        self.width = 100
        self.length = 100*4*10
        self.gate_width = 100
        self.gate_length = 300
        self.n_clav_rep = 10
        self.spacing = 100*4
        self.shift = 0
        self.x = 0
        self.y = 0
        self.fillet = 0  # 0.02
        self.fillet_tolerance = 1e-4
        self.rotation = 0

    def build(self):
        """
        Construct the Clavier gate's geometric representation based on the specified properties.

        This method generates the detailed geometry of the Clavier gate, taking into account its width, length, number of repetitions, and other characteristics. The process involves creating a complex polygon structure, performing necessary geometric transformations, and finalizing the gate's shape.

        The constructed Clavier gate is then registered within a cell, integrating its geometry into the overall quantum dot device design.

        Raises:
            ValueError: If no valid Layer is assigned before building.
        """
        cl_points = generate_clavier_gates(self.width, self.length,
                                           self.gate_width,
                                           self.gate_length,
                                           self.n_clav_rep, self.spacing,
                                           self.shift, (-self.length/2, 0),
                                           self.rotation)
        layer = getattr(self.layer, self.layer_stage)
        cl = gdstk.Polygon(cl_points, layer=layer)

        cl.translate(self.x, self.y)
        cl.fillet(self.fillet, tolerance=self.fillet_tolerance)
        cell = gdstk.Cell(self.name)
        cell.add(cl)
        self.elements[self.name]['vertices'] = cl.points
        self.elements[self.name]['positions'] = [[self.x, self.y]]
        self.elements[self.name]['layer'] = self.layer
        self.elements[self.name]['layer_stage'] = self.layer_stage
        self.cell = cell
        self._set_built(True)
