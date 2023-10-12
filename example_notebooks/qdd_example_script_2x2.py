# -*- coding: utf-8 -*-
"""
Created on Tue May  9 08:18:02 2023

@author: vjohn
"""
# %% import

import numpy as np
import QuantumDotDesigner as qdd

from QuantumDotDesigner.base import Layer

from QuantumDotDesigner.elements import (Plunger, Barrier, ScreeningGate,
                                         BasicPolygon)
from QuantumDotDesigner.components import Sensor, FanOutLine

from QuantumDotDesigner.helpers.helpers import mirror_points_along_path as mirror

# %% Elements/components

collection = qdd.BaseCollection()

# %% Layers

ohmic_layer = Layer('ohmic_layer', collection)
ohmic_layer.fine = 3
ohmic_layer.coarse = 4

barrier_layer = Layer('barrier_layer', collection)
barrier_layer.fine = 5
barrier_layer.coarse = 6

screening_layer = Layer('screening_layer', collection, 3, 4)
screening_layer.fine = 31
screening_layer.coarse = 32

barrier_source_layer = barrier_layer.copy('barrier_source_layer', collection)
barrier_drain_layer = barrier_layer.copy('barrier_drain_layer', collection)
barrier_sep_layer = screening_layer.copy('barrier_screening_layer', collection)

plunger_layer = Layer('plunger_layer', collection)
plunger_layer.fine = 21
plunger_layer.coarse = 22
plunger_layer.via_etch = 23
plunger_layer.via_fine = 24
plunger_layer.via_coarse = 25

# %%% define plungers

pl = Plunger('plunger', collection)
pl.diameter = 150e-3
pl.layer = plunger_layer

# %%% define barriers

bar_0deg = Barrier('barrier_0deg_rotated', collection)
bar_0deg.width = 40e-3
bar_0deg.length = 70e-3
bar_0deg.layer = barrier_layer
bar_0deg.rotate = 0/2*np.pi

bar_90deg = bar_0deg.copy('barrier_90deg_rotated', collection)
bar_90deg.rotate = 1/2*np.pi

bar_180deg = bar_0deg.copy('barrier_180deg_rotated', collection)
bar_180deg.rotate = 2/2*np.pi

bar_270deg = bar_0deg.copy('barrier_270deg_rotated', collection)
bar_270deg.rotate = 3/2*np.pi

# %%% define sensor

spacing_sep = 60e-3

sensor_top = Sensor('sensor_top', collection)

# for the positioning of the sensor element you can either indicate 'top',
# 'bottom', 'top-right', etc ...., or you indiacte the angle with 'top'
# corresponding to 0, 'top-right' corresponding to np.pi/4, etc ...
sensor_top.sep_pos = 'bottom-right'
# sensor_top.sep_pos_angle = 1/8*np.pi
sensor_top.source_pos = 'top-right'
# sensor_top.source_pos_angle = 3/8*np.pi
sensor_top.drain_pos = 'bottom-left'  # 1/8*np.pi
# sensor_top.drain_pos_angle = -3/8*np.pi

# the orientation of the sensor elements is indicated by the direction in
# which the end of the element points to, This can be either clockwise or
# counterclockwise
sensor_top.barrier_orientation['drain'] = 'counterclockwise'
sensor_top.barrier_orientation['source'] = 'clockwise'
sensor_top.barrier_orientation['sep'] = 'clockwise'

sensor_top.gap_sep = 60e-3
sensor_top.gap_ohmic_pl = 50e-3

sensor_top.plunger.diameter = 160e-3
sensor_top.plunger.layer = plunger_layer

sensor_top.barrier_source.width = 40e-3
sensor_top.barrier_source.length = 70e-3
sensor_top.barrier_source.layer = barrier_source_layer

sensor_top.barrier_drain.width = 40e-3
sensor_top.barrier_drain.length = 70e-3
sensor_top.barrier_drain.layer = barrier_drain_layer

sensor_top.barrier_sep.width = 50e-3
sensor_top.barrier_sep.length = 60e-3
sensor_top.barrier_sep.layer = screening_layer

sensor_top.source.layer = ohmic_layer
sensor_top.source.contact_length = 70e-3
sensor_top.drain.layer = ohmic_layer
sensor_top.drain.contact_length = 70e-3

sensor_bottom = sensor_top.copy('sensor_bottom', collection)
sensor_bottom.sep_pos = 'top-left'
sensor_bottom.source_pos = 'bottom-left'
sensor_bottom.drain_pos = 'top-right'

# sensor_top.plot(build=True)

# %% Unit cell

unit_cell = qdd.UnitCell('unit_cell')
spacing_qd = 200e-3

# %%% add plunger
uc_pl = unit_cell.add_component()

uc_bar_0deg = unit_cell.add_component()
uc_bar_90deg = unit_cell.add_component()
uc_bar_180deg = unit_cell.add_component()
uc_bar_270deg = unit_cell.add_component()

uc_pl.component = pl
uc_pl.center = (0, 0)
uc_pl.rows = 2
uc_pl.columns = 2
uc_pl.spacing = (spacing_qd, spacing_qd)

# %%% add barriers

uc_bar_0deg.component = bar_0deg
uc_bar_0deg.center = (spacing_qd/2, 0)

uc_bar_90deg.component = bar_90deg
uc_bar_90deg.center = (0, spacing_qd/2)

uc_bar_180deg.component = bar_180deg
uc_bar_180deg.center = (-spacing_qd/2, 0)

uc_bar_270deg.component = bar_270deg
uc_bar_270deg.center = (0, -spacing_qd/2)

# %%% add sensor

sensor_pos_uniy = (spacing_qd/2 +
                   2**0.5/2*(pl.diameter/2 * pl.asymy +
                             sensor_top.gap_sep +
                             sensor_top.plunger.diameter/2 *
                             sensor_top.plunger.asymy))

sensor_pos_unix = (spacing_qd/2 +
                   2**0.5/2*(pl.diameter/2 * pl.asymx +
                             sensor_top.gap_sep +
                             sensor_top.plunger.diameter/2 *
                             sensor_top.plunger.asymx))

uc_st = unit_cell.add_component()
uc_st.component = sensor_top
uc_st.center = (-sensor_pos_unix, sensor_pos_uniy)

uc_sb = unit_cell.add_component()
uc_sb.component = sensor_bottom
uc_sb.center = (sensor_pos_unix, -sensor_pos_uniy)


unit_cell.build()
unit_cell.plot()

# %% Main cell

qda = qdd.QuantumDotArray()

# %%% add unit cell

uc_unitcell = qda.add_component()
uc_unitcell.component = unit_cell

# %% Fanout

fo_points = qdd.FanoutPoints(qda)
fo = qdd.Fanout('fanout')

fo_points.fanout_counts = {'top': 6, 'bottom': 6, 'left': 7, 'right': 7}

fo_points.fo_stages = [(16, 16), (500, 530), (1200, 1200)]
# fo_points.bondpad_position = {'top':  1500, 'bottom': 1500,
#                               'left': 1500, 'right': 1500}
fo_points.fo_widths = [1, 6, 25]
fo_points.spacings = [2, 80, 250]

fo_points.bondpad_size = {'top': (110, 400), 'bottom': (110, 400),
                          'left': (400, 110), 'right': (400, 110)}

fo_points.ohmic_bondpad['top']['shift'] = 50
fo_points.ohmic_bondpad['bottom']['shift'] = 50
fo_points.ohmic_bondpad['left']['shift'] = 90
fo_points.ohmic_bondpad['right']['shift'] = 90

fo_points.create_fo_polygons_coarse()

# %%% Ohmic fanout parameters

bp_ohmic_position = 1000
bp_ohmic_width_out = 110
bp_ohmic_width_in = 80
bp_ohmic_length = 400

bp_shift_top_source = 50
bp_shift_top_drain = 100
bp_shift_bottom_source = 50
bp_shift_bottom_drain = 100

# %%% Fanout plunger
# %%%% Fanout plunger ver 0
fo_pl_0 = FanOutLine('plunger', 0, collection, fo_points)

fo_pl_0.fo_direction = 'left'
fo_pl_0.n_fanout = 4

fo_pl_0.fo_line_fine.fo_width_start = 40e-3
fo_pl_0.fo_line_fine.points_along_path = [[-0.05, -0.01, 'start'],
                                          [-0.4, -0.1, 'prev']
                                          ]

fo.add_component(fo_pl_0)

# %%%% Fanout plunger ver 1
fo_pl_1 = FanOutLine('plunger', 1, collection, fo_points)

fo_pl_1.fo_direction = 'right'
fo_pl_1.n_fanout = 0

fo_pl_1.fo_line_fine.fo_width_start = 40e-3
fo_pl_1.fo_line_fine.points_along_path = [[0.6, 0.6, 'start']
                                          ]

fo.add_component(fo_pl_1)

# %%%% Fanout plunger ver 2
fo_pl_2 = FanOutLine('plunger', 2, collection, fo_points)

fo_pl_2.fo_direction = 'left'
fo_pl_2.n_fanout = -1

fo_pl_2.fo_line_fine.fo_width_start = 40e-3
points = mirror(fo_pl_1.fo_line_fine.points_along_path)
fo_pl_2.fo_line_fine.points_along_path = points

fo.add_component(fo_pl_2)

# %%%% Fanout plunger ver 3
fo_pl_3 = FanOutLine('plunger', 3, collection, fo_points)

fo_pl_3.fo_direction = 'right'
fo_pl_3.n_fanout = 2

fo_pl_3.fo_line_fine.fo_width_start = 40e-3
points = mirror(fo_pl_0.fo_line_fine.points_along_path)
fo_pl_3.fo_line_fine.points_along_path = points

fo.add_component(fo_pl_3)

# %%% Fanout barriers
# %%%% Fanout bar 0deg 0
fo_bar_0deg_0 = FanOutLine('barrier_0deg_rotated', 0, collection, fo_points)

fo_bar_0deg_0.fo_direction = 'right'
fo_bar_0deg_0.n_fanout = 1

fo_bar_0deg_0.fo_line_fine.fo_width_start = 40e-3
fo_bar_0deg_0.fo_line_fine.points_along_path = [[0.1, 0, 'start'],
                                                [0.1, 0.03, 'prev'],
                                                [0.2, 0.1, 'prev']
                                                ]

fo.add_component(fo_bar_0deg_0)

# %%%% Fanout bar 90deg 0
fo_bar_90deg_0 = FanOutLine('barrier_90deg_rotated', 0, collection, fo_points)

fo_bar_90deg_0.fo_direction = 'top'
fo_bar_90deg_0.n_fanout = 4

fo_bar_90deg_0.fo_line_fine.fo_width_start = 40e-3
fo_bar_90deg_0.fo_line_fine.points_along_path = [[0, 0.1, 'start'],
                                                 [0.03, 0.1, 'prev'],
                                                 [0.1, 0.2, 'prev']
                                                 ]

fo.add_component(fo_bar_90deg_0)

# %%%% Fanout bar 180deg 0
fo_bar_180deg_0 = FanOutLine(
    'barrier_180deg_rotated', 0, collection, fo_points)

fo_bar_180deg_0.fo_direction = 'left'
fo_bar_180deg_0.n_fanout = 5

fo_bar_180deg_0.fo_line_fine.fo_width_start = 40e-3
points = mirror(fo_bar_0deg_0.fo_line_fine.points_along_path)
fo_bar_180deg_0.fo_line_fine.points_along_path = points

fo.add_component(fo_bar_180deg_0)

# %%%% Fanout bar 270deg 0
fo_bar_270deg_0 = FanOutLine(
    'barrier_270deg_rotated', 0, collection, fo_points)

fo_bar_270deg_0.fo_direction = 'bottom'
fo_bar_270deg_0.n_fanout = 1

fo_bar_270deg_0.fo_line_fine.fo_width_start = 40e-3
points = mirror(fo_bar_90deg_0.fo_line_fine.points_along_path)
fo_bar_270deg_0.fo_line_fine.points_along_path = points

fo.add_component(fo_bar_270deg_0)

# %%% Fanout sensor top
# %%%% Fanout sensor plunger top
fo_sens_pl_top = FanOutLine('sensor_top_plunger', 0, collection, fo_points)

fo_sens_pl_top.fo_direction = 'left'
fo_sens_pl_top.n_fanout = 0

fo_sens_pl_top.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_top.fo_line_fine.points_along_path = [[-0.5, 0.5, 'start', 40e-3]
                                                 ]

fo.add_component(fo_sens_pl_top)

# %%%% Fanout sensor source top
fo_sens_top_source = FanOutLine(
    'sensor_top_source', 0, collection, fo_points)

fo_sens_top_source.fo_direction = 'top'
fo_sens_top_source.n_fanout = 2

fo_sens_top_source.bp_ohmic_position = bp_ohmic_position
fo_sens_top_source.bp_ohmic_width_out = bp_ohmic_width_out
fo_sens_top_source.bp_ohmic_width_in = bp_ohmic_width_in
fo_sens_top_source.bp_ohmic_length = bp_ohmic_length
fo_sens_top_source.bp_ohmic_shift = bp_shift_top_source

fo_sens_top_source.fo_line_fine.points_along_path = [[0, 0.4, 'start'],
                                                     [-0.04, 0.4, 'prev']
                                                     ]

fo.add_component(fo_sens_top_source)

# %%%% Fanout sensor drain top
fo_sens_top_drain = FanOutLine(
    'sensor_top_drain', 0, collection, fo_points)

fo_sens_top_drain.fo_direction = 'left'
fo_sens_top_drain.n_fanout = 2

fo_sens_top_drain.bp_ohmic_position = bp_ohmic_position
fo_sens_top_drain.bp_ohmic_width_out = bp_ohmic_width_out
fo_sens_top_drain.bp_ohmic_width_in = bp_ohmic_width_in
fo_sens_top_drain.bp_ohmic_length = bp_ohmic_length
fo_sens_top_drain.bp_ohmic_shift = bp_shift_top_drain

fo_sens_top_drain.fo_line_fine.points_along_path = [[-0.25, -0.01, 'start']
                                                    ]

fo.add_component(fo_sens_top_drain)

# %%%% Fanout sensor source barrier top
fo_sens_top_bar_source = FanOutLine(
    'sensor_top_barrier_source', 0, collection, fo_points)

fo_sens_top_bar_source.fo_direction = 'top'
fo_sens_top_bar_source.n_fanout = 1

fo_sens_top_bar_source.fo_line_fine.points_along_path = [[-0.1, 0.1, 'start'],
                                                         [-0.1, 0.2, 'prev']
                                                         ]

fo.add_component(fo_sens_top_bar_source)

# %%%% Fanout sensor drain barrier top
fo_sens_top_bar_drain = FanOutLine(
    'sensor_top_barrier_drain', 0, collection, fo_points)

fo_sens_top_bar_drain.fo_direction = 'left'
fo_sens_top_bar_drain.n_fanout = 1

fo_sens_top_bar_drain.fo_line_fine.points_along_path = [[-0.1, 0.1, 'start'],
                                                        [-0.2, 0.1, 'prev']
                                                        ]

fo.add_component(fo_sens_top_bar_drain)

# %%%% Fanout sensor sep barrier top
fo_sens_top_bar_sep = FanOutLine(
    'sensor_top_barrier_seperation', 0, collection, fo_points)

fo_sens_top_bar_sep.fo_direction = 'top'
fo_sens_top_bar_sep.n_fanout = 3

fo_sens_top_bar_sep.fo_line_fine.points_along_path = [[0.1, 0.1, 'start'],
                                                      [0.1, 0.2, 'prev'],
                                                      [0, 0.3, 'prev']
                                                      ]

fo.add_component(fo_sens_top_bar_sep)

# %%% Fanout sensor bottom
# %%%% Fanout sensor plunger bottom
fo_sens_pl_bottom = FanOutLine(
    'sensor_bottom_plunger', 0, collection, fo_points)

fo_sens_pl_bottom.fo_direction = 'right'
fo_sens_pl_bottom.n_fanout = -1

fo_sens_pl_bottom.fo_line_fine.fo_width_start = 40e-3
points = mirror(fo_sens_pl_top.fo_line_fine.points_along_path)
fo_sens_pl_bottom.fo_line_fine.points_along_path = points

fo.add_component(fo_sens_pl_bottom)

# %%%% Fanout sensor source bottom
fo_sens_bottom_source = FanOutLine(
    'sensor_bottom_source', 0, collection, fo_points)

fo_sens_bottom_source.fo_direction = 'bottom'
fo_sens_bottom_source.n_fanout = 3

fo_sens_bottom_source.bp_ohmic_position = bp_ohmic_position
fo_sens_bottom_source.bp_ohmic_width_out = bp_ohmic_width_out
fo_sens_bottom_source.bp_ohmic_width_in = bp_ohmic_width_in
fo_sens_bottom_source.bp_ohmic_length = bp_ohmic_length
fo_sens_bottom_source.bp_ohmic_shift = bp_shift_bottom_source

points = mirror(fo_sens_top_source.fo_line_fine.points_along_path)
fo_sens_bottom_source.fo_line_fine.points_along_path = points

fo.add_component(fo_sens_bottom_source)

# %%%% Fanout sensor drain bottom
fo_sens_bottom_drain = FanOutLine(
    'sensor_bottom_drain', 0, collection, fo_points)

fo_sens_bottom_drain.fo_direction = 'right'
fo_sens_bottom_drain.n_fanout = 4

fo_sens_bottom_drain.bp_ohmic_position = bp_ohmic_position
fo_sens_bottom_drain.bp_ohmic_width_out = bp_ohmic_width_out
fo_sens_bottom_drain.bp_ohmic_width_in = bp_ohmic_width_in
fo_sens_bottom_drain.bp_ohmic_length = bp_ohmic_length
fo_sens_bottom_drain.bp_ohmic_shift = bp_shift_bottom_drain

points = mirror(fo_sens_top_drain.fo_line_fine.points_along_path)
fo_sens_bottom_drain.fo_line_fine.points_along_path = points

fo.add_component(fo_sens_bottom_drain)

# %%%% Fanout sensor source barrier bottom
fo_sens_bottom_bar_source = FanOutLine(
    'sensor_bottom_barrier_source', 0, collection, fo_points)

fo_sens_bottom_bar_source.fo_direction = 'bottom'
fo_sens_bottom_bar_source.n_fanout = 4

points = mirror(fo_sens_top_bar_source.fo_line_fine.points_along_path)
fo_sens_bottom_bar_source.fo_line_fine.points_along_path = points

fo.add_component(fo_sens_bottom_bar_source)

# %%%% Fanout sensor drain barrier bottom
fo_sens_bottom_bar_drain = FanOutLine(
    'sensor_bottom_barrier_drain', 0, collection, fo_points)

fo_sens_bottom_bar_drain.fo_direction = 'right'
fo_sens_bottom_bar_drain.n_fanout = 5

points = mirror(fo_sens_top_bar_drain.fo_line_fine.points_along_path)
fo_sens_bottom_bar_drain.fo_line_fine.points_along_path = points

fo.add_component(fo_sens_bottom_bar_drain)

# %%%% Fanout sensor sep barrier bottom
fo_sens_bottom_bar_sep = FanOutLine(
    'sensor_bottom_barrier_seperation', 0, collection, fo_points)

fo_sens_bottom_bar_sep.fo_direction = 'bottom'
fo_sens_bottom_bar_sep.n_fanout = 2

points = mirror(fo_sens_top_bar_sep.fo_line_fine.points_along_path)
fo_sens_bottom_bar_sep.fo_line_fine.points_along_path = points

fo.add_component(fo_sens_bottom_bar_sep)

# %% Screening gates
# %%% Plunger screening
# %%%% Plunger 0 screening
screen_pl_0 = ScreeningGate('screening_gate_pl_0', collection)

screen_pl_0.screen('plunger', 0,
                   [0.1, 0.3, 0.4], [50e-3, 50e-3, (25e-3, 75e-3)])
screen_pl_0.layer = screening_layer
screen_pl_0.fo_contact_width = 50e-3
screen_pl_0.fo_contact_direction = 1

screen_pl_0_qda = qda.add_component()
screen_pl_0_qda.component = screen_pl_0

# %%%% Plunger 1 screening
screen_pl_1 = ScreeningGate('screening_gate_pl_1', collection)

screen_pl_1.screen('plunger', 1,
                   [0.1, 0.3, 0.4], [50e-3, 50e-3, (75e-3, 25e-3)])
screen_pl_1.layer = screening_layer
screen_pl_1.fo_contact_width = 50e-3
screen_pl_1.fo_contact_direction = 0

screen_pl_1_qda = qda.add_component()
screen_pl_1_qda.component = screen_pl_1

# %%%% Plunger 2 screening
screen_pl_2 = ScreeningGate('screening_gate_pl_2', collection)

screen_pl_2.screen('plunger', 2,
                   [0.1, 0.3, 0.4], [50e-3, 50e-3, (75e-3, 25e-3)])
screen_pl_2.layer = screening_layer
screen_pl_2.fo_contact_width = 50e-3
screen_pl_2.fo_contact_direction = 0

screen_pl_2_qda = qda.add_component()
screen_pl_2_qda.component = screen_pl_2

# %%%% Plunger 3 screening
screen_pl_3 = ScreeningGate('screening_gate_pl_3', collection)

screen_pl_3.screen('plunger', 3,
                   [0.1, 0.3, 0.4], [50e-3, 50e-3, (25e-3, 75e-3)])
screen_pl_3.layer = screening_layer
screen_pl_3.fo_contact_width = 50e-3
screen_pl_3.fo_contact_direction = 1

screen_pl_3_qda = qda.add_component()
screen_pl_3_qda.component = screen_pl_3

# %%% Sensor plunger screening
# %%%% Plunger sens top plunger screening
screen_sens_pl_top = ScreeningGate('screening_gate_sens_pl_top',
                                   collection)

screen_sens_pl_top.screen('sensor_top_plunger', 0,
                          [0.1, 0.3, 0.4], [50e-3, 50e-3, (25e-3, 75e-3)])
screen_sens_pl_top.layer = screening_layer
screen_sens_pl_top.fo_contact_width = 50e-3
screen_sens_pl_top.fo_contact_direction = 1

screen_sens_pl_top_qda = qda.add_component()
screen_sens_pl_top_qda.component = screen_sens_pl_top

# %%%% Plunger sens bottom plunger screening
screen_sens_pl_bottom = ScreeningGate('screening_gate_sens_pl_bottom',
                                      collection)

screen_sens_pl_bottom.screen('sensor_bottom_plunger', 0,
                             [0.1, 0.3, 0.4], [50e-3, 50e-3, (25e-3, 75e-3)])
screen_sens_pl_bottom.layer = screening_layer
screen_sens_pl_bottom.fo_contact_width = 50e-3
screen_sens_pl_bottom.fo_contact_direction = 1

screen_sens_pl_bottom_qda = qda.add_component()
screen_sens_pl_bottom_qda.component = screen_sens_pl_bottom

# %% Fanout screening gates
qda.build()

# %%% Plunger screening fanout
# %%%% Fanout screening plunger 0

fo_screen_pl_0 = FanOutLine('screening_gate_pl_0', 0, collection, fo_points)

fo_screen_pl_0.fo_direction = 'left'
fo_screen_pl_0.n_fanout = 3

fo_screen_pl_0.fo_line_fine.fo_width_start = screen_pl_0.fo_contact_width
fo_screen_pl_0.fo_line_fine.points_along_path = [[-0.1, 0, 'start'],
                                                 [-0.8, 0.01, 'start'],
                                                 ]

fo.add_component(fo_screen_pl_0)

# %%%% Fanout screening plunger 1

fo_screen_pl_1 = FanOutLine(
    'screening_gate_pl_1', 0, collection, fo_points)

fo_screen_pl_1.fo_direction = 'top'
fo_screen_pl_1.n_fanout = -1

fo_screen_pl_1.fo_line_fine.fo_width_start = screen_pl_1.fo_contact_width
fo_screen_pl_1.fo_line_fine.points_along_path = [[0.05, 0.1, 'start'],
                                                 [0.1, 0.3, 'start'],
                                                 ]

fo.add_component(fo_screen_pl_1)

# %%%% Fanout screening plunger 2

fo_screen_pl_2 = FanOutLine(
    'screening_gate_pl_2', 0, collection, fo_points)

fo_screen_pl_2.fo_direction = 'bottom'
fo_screen_pl_2.n_fanout = 0

fo_screen_pl_2.fo_line_fine.fo_width_start = screen_pl_2.fo_contact_width
points = mirror(fo_screen_pl_1.fo_line_fine.points_along_path)
fo_screen_pl_2.fo_line_fine.points_along_path = points

fo.add_component(fo_screen_pl_2)

# %%%% Fanout screening plunger 3

fo_screen_pl_3 = FanOutLine(
    'screening_gate_pl_3', 0, collection, fo_points)

fo_screen_pl_3.fo_direction = 'right'
fo_screen_pl_3.n_fanout = 3

fo_screen_pl_3.fo_line_fine.fo_width_start = screen_pl_3.fo_contact_width
points = mirror(fo_screen_pl_0.fo_line_fine.points_along_path)
fo_screen_pl_3.fo_line_fine.points_along_path = points

fo.add_component(fo_screen_pl_3)

# %%% Sensor plunger screening fanout
# %%%% Fanout screening plunger sensor top

fo_screen_sens_pl_top = FanOutLine(
    'screening_gate_sens_pl_top', 0, collection, fo_points)

fo_screen_sens_pl_top.fo_direction = 'top'
fo_screen_sens_pl_top.n_fanout = 0

fo_screen_sens_pl_top.fo_line_fine.fo_width_start = screen_sens_pl_top.fo_contact_width
fo_screen_sens_pl_top.fo_line_fine.points_along_path = [[-0.05, 0.1, 'start'],
                                                        [-0.2, 0.4, 'start'],
                                                        ]

fo.add_component(fo_screen_sens_pl_top)

# %%%% Fanout screening plunger sensor bottom

fo_screen_sens_pl_bottom = FanOutLine(
    'screening_gate_sens_pl_bottom', 0, collection, fo_points)

fo_screen_sens_pl_bottom.fo_direction = 'bottom'
fo_screen_sens_pl_bottom.n_fanout = -1

fo_screen_sens_pl_bottom.fo_line_fine.fo_width_start = screen_sens_pl_bottom.fo_contact_width
points = mirror(fo_screen_sens_pl_top.fo_line_fine.points_along_path)
fo_screen_sens_pl_bottom.fo_line_fine.points_along_path = points

fo.add_component(fo_screen_sens_pl_bottom)

# %% Main cell with fanout

fo_qda = qda.add_component()
fo_qda.component = fo

# %% Build and save

qda.build()
layout_path = "..\example_notebooks\layout_files\chip_layout.gds"
layout = qda.add_chip_layout(layout_path)
qda.save_as_gds('..\example_notebooks\example_2x2_device.gds')
