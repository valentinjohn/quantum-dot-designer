# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:02:28 2023

@author: vjohn
"""

from QuantumDotDesigner.elements.FanOutLineBase import FanOutLineBase


class FanOutLineCoarse(FanOutLineBase):
    pass


# from QuantumDotDesigner.elements.FanOutLineBase import FanOutLineBase
# from QuantumDotDesigner.base import BaseCollection
# import gdstk


# class FanOutLineCoarse(FanOutLineBase):
#     def __init__(self, name, collection: BaseCollection):
#         super().__init__(name, collection)

#     def build(self):
#         fo_line = gdstk.Polygon(self.polygons, layer=self.layer.coarse)
#         fo_line.fillet(self.fillet, tolerance=self.fillet_tolerance)

#         self.elements[self.name]['vertices'] = self.polygons  # fo_line.points
#         self.elements[self.name]['positions'] = [[0, 0]]
#         self.elements[self.name]['layer'] = self.layer

#         self.cell.add(fo_line)
#         self._set_built(True)
