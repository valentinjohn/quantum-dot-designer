# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 10:46:45 2023

@author: vjohn
"""

import gdstk
from QuantumDotDesigner.helpers.helpers import merge_device_positions
from QuantumDotDesigner.base import UnitCell
from QuantumDotDesigner.base import PlotMixin, Sublattice


class QuantumDotArray(PlotMixin):
    def __init__(self):
        """
        Initialize a QuantumDotArray object.

        Args:
            parent_instance: The parent instance.
            unitcell_instance (UnitCell): An optional UnitCell instance to use
            for the array.
        """
        self.name = 'Quantum_Dot_Array'
        self.spacing_qd = 200e-3
        self.spacing_qd_diag = 2**0.5 * self.spacing_qd
        self.elements = {}
        self.components = {}
        self.components_position = {}
        self.main_cell = gdstk.Cell('MAIN')
        self._n = 0
        self.chip_layout_path = None
        self._built = False

    @property
    def built(self):
        return self._built

    def add_component(self):
        """
        Add a sublattice to the QuantumDotArray.

        Args:
            name (str): Name of the sublattice.

        Returns:
            Sublattice: The created Sublattice object.
        """
        name = f'MAIN_sublattice_{self._n}'
        self._n = self._n + 1
        sublattice = Sublattice(name)
        self.components[name] = sublattice
        return sublattice

    def add_chip_layout(self, chip_layout_path=None):
        if chip_layout_path is not None:
            self.chip_layout_path = chip_layout_path
        layout = gdstk.read_rawcells(chip_layout_path)

        self.main_cell.add(gdstk.Reference(layout['TOP']))
        self.main_cell.flatten()

        return layout

    def build(self):
        """
        Build the QuantumDotArray.

        This method adds the sublattices and unit cells to the main cell object.
        """
        elements = {}
        for cell in self.components.values():
            if not cell.built:
                cell.build()
            elements = merge_device_positions(elements, cell.elements)
            if isinstance(cell, Sublattice):
                self.main_cell.add(gdstk.Reference(cell.cell))
            elif isinstance(cell, UnitCell):
                for c in cell.cells:
                    self.main_cell.add(gdstk.Reference(c.cell))
        self.main_cell.flatten()
        self.elements = elements
        self._built = True

    def save_as_gds(self, filename):
        """
        Save the QuantumDotArray as a GDS file.

        Args:
            filename (str): Name of the output GDS file.
        """
        lib = gdstk.Library()
        lib.add(self.main_cell, *self.main_cell.dependencies(True))
        lib.write_gds(filename)
