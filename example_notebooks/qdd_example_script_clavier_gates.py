# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 16:28:10 2023

@author: vjohn
"""

# %% import

import numpy as np
import quantum_dot_designer as qdd

from quantum_dot_designer.base.Layer import Layer
from quantum_dot_designer.elements import Plunger, Barrier, ScreeningGate
from quantum_dot_designer.components import Sensor, FanOutLine, Clavier

# %% Init

collection = qdd.BaseCollection()
unit_cell = qdd.UnitCell('unit_cell')
qda = qdd.QuantumDotArray()

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

clav_gate_1_layer = barrier_layer.copy('clav_gate_1_layer', collection)
clav_gate_2_layer = screening_layer.copy('clav_gate_2_layer', collection)

# %% define plungers

pl = Plunger('plunger', collection)
pl.asymx = 1
pl.diameter = 130e-3
pl.layer = plunger_layer

# pl.plot()

# %% define barriers

bars_rot = np.pi*(-1/2+np.array([0, -2/3, 2/3]))

bar_1 = Barrier('barrier_1', collection)
bar_1.width = 40e-3
bar_1.length = 30e-3
bar_1.layer = barrier_layer
bar_1.rotate = bars_rot[0]

bar_2 = bar_1.copy('barrier_2', collection)
bar_2.rotate = bars_rot[1]

bar_3 = bar_1.copy('barrier_3', collection)
bar_3.rotate = bars_rot[2]

# bar_1.plot()
# bar_2.plot()
# bar_3.plot()

# %% define unit cell

spacing_qd = 200e-3
tritop_height = 3**0.5/2 * spacing_qd

uc_pl_tribase = unit_cell.add_component()
uc_pl_tritop = unit_cell.add_component()

uc_bar_1 = unit_cell.add_component()
uc_bar_2 = unit_cell.add_component()
uc_bar_3 = unit_cell.add_component()

uc_pl_tribase.component = pl
uc_pl_tribase.center = (0, 0)
uc_pl_tribase.columns = 2
uc_pl_tribase.spacing = (qda.spacing_qd,
                         qda.spacing_qd)

uc_pl_tritop.component = pl
uc_pl_tritop.center = (0, tritop_height)

uc_bar_1.component = bar_1
uc_bar_1.center = (0, 0)

uc_bar_2.component = bar_2
uc_bar_2.center = (-1*spacing_qd/4, tritop_height/2)

uc_bar_3.component = bar_3
uc_bar_3.center = (1*spacing_qd/4, tritop_height/2)

# %% define sensor

sensor_top = Sensor('sensor_top', collection)

sensor_top.sep_pos = 'bottom'
sensor_top.source_pos = 'left'
sensor_top.drain_pos = 'right'

sensor_top.bar_drain_end = 'clockwise'
sensor_top.bar_source_end = 'counterclockwise'
sensor_top.bar_sep_end = 'clockwise'

sensor_top.gap_sep = 60e-3
sensor_top.gap_ohmic_pl = 40e-3

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
sensor_top.barrier_sep.layer = barrier_layer

sensor_top.source.layer = ohmic_layer
sensor_top.drain.layer = ohmic_layer

# sensor_top.plot()

# %% add sensor to unit cell

unit_cell_ylim = (tritop_height +
                  pl.diameter/2 *
                  pl.asymy)

sensor_pos_uniy = (unit_cell_ylim+sensor_top.gap_sep +
                   sensor_top.plunger.diameter/2 *
                   sensor_top.plunger.asymy)

sensor_pos_unix = 0

uc_st = unit_cell.add_component()
uc_st.component = sensor_top
uc_st.center = (sensor_pos_unix, sensor_pos_uniy)

# %% Add unit cell to main cell

# unit_cell.plot()

uc_unitcell = qda.add_component()
uc_unitcell.component = unit_cell
uc_unitcell.center = (0, 0)
uc_unitcell.rows = 1
uc_unitcell.columns = 1
uc_unitcell.spacing = (2*qda.spacing_qd_diag,
                       qda.spacing_qd_diag)
# %% Add clavier gates

clavier = Clavier('clavier', collection)

clavier.clav_dot_size = 100e-3
clavier.clav_gate_gap = 20e-3
clavier.clav_gap = [.4, .7]
clavier.clav_width = 200e-3
clavier.n_clav_rep = 8

clavier.screen_width = 100e-3
clavier.screen_gap = 0
clavier.screen_layer = screening_layer

clavier.clav_layers = [clav_gate_1_layer, clav_gate_2_layer]
clavier.x = -(spacing_qd + clavier.screen_length/2)
clavier.y = 0
clavier.fillet = 0.02
clavier.fillet_tolerance = 1e-4

clavier.mirror = False

clavier.build()

clavier_mirrored = clavier.copy('clavier_mirrored', collection)
clavier_mirrored.mirror = True
clavier_mirrored.x = (spacing_qd + clavier.screen_length/2)

clavier_mirrored.build()

sl_clavier = qda.add_component()
sl_clavier.component = clavier

sl_clavier_mirrored = qda.add_component()
sl_clavier_mirrored.component = clavier_mirrored

# clavier.plot(build=True)

# %% Fanout

fo_points = qdd.FanoutPoints(qda)
fo = qdd.Fanout('fanout')

fo_points.fanout_counts = {'top': 10, 'bottom': 10, 'left': 5, 'right': 5}

fo_points.fo_stages = [(16, 16), (500, 530), (1200, 1200)]
fo_points.bondpad_position = {'top':  1500, 'bottom': 1500,
                              'left': 1500, 'right': 1500}
fo_points.fo_widths = [1, 6, 25]
fo_points.spacings = [2, 80, 250]

fo_points.bondpad_size = {'top': (110, 400), 'bottom': (110, 400),
                          'left': (400, 110), 'right': (400, 110)}
fo_points.create_fo_polygons_coarse()

# %%% Fanout plungers
# %%%% Fanout plunger 0

fo_pl_0 = FanOutLine('plunger', 0, collection, fo_points)

fo_pl_0.fo_direction = 'top'
fo_pl_0.n_fanout = 8

fo_pl_0.fo_line_fine.fine_fo_width_start = 40e-3
fo_pl_0.fo_line_fine.points_along_path = [[0.25, 0.15, 'start'],
                                          [0.15, 0.8, 'prev']]

fo.add_component(fo_pl_0)

# %%%% Fanout plunger 1
fo_pl_1 = FanOutLine('plunger', 1, collection, fo_points)

fo_pl_1.fo_direction = 'bottom'
fo_pl_1.n_fanout = 4

fo_pl_1.fo_line_fine.fo_width_start = 40e-3
fo_pl_1.fo_line_fine.points_along_path = [[0, -0.1, 'start'],
                                          [0, -0.7, 'prev']]
fo.add_component(fo_pl_1)

# %%%% Fanout plunger 2
fo_pl_2 = FanOutLine('plunger', 2, collection, fo_points)

fo_pl_2.fo_direction = 'bottom'
fo_pl_2.n_fanout = 6

fo_pl_2.fo_line_fine.fo_width_start = 40e-3
fo_pl_2.fo_line_fine.points_along_path = [[0, -0.1, 'start'],
                                          [0, -0.7, 'prev']]
fo.add_component(fo_pl_2)

# %%% Fanout barriers
# %%%% Fanout bar 1
fo_bar_1 = FanOutLine('barrier_1', 0, collection, fo_points)


fo_bar_1.fo_direction = 'bottom'
fo_bar_1.n_fanout = 5

fo_bar_1.fo_line_fine.fo_width_start = 40e-3
fo_bar_1.fo_line_fine.points_along_path = [[0, -0.5, 'start'],
                                           # [-0.1, 0.1, 'prev'],
                                           # [0, 0.3, 'prev'],
                                           # [-0.15, 0.3, 'prev']
                                           ]
fo.add_component(fo_bar_1)

# %%%% Fanout bar 2
fo_bar_2 = FanOutLine('barrier_2', 0, collection, fo_points)

fo_bar_2.fo_direction = 'top'
fo_bar_2.n_fanout = 1

fo_bar_2.fo_line_fine.fo_width_start = 40e-3
fo_bar_2.fo_line_fine.points_along_path = [[-0.1, 0.1, 'start'],
                                           [-0.15, 0.05, 'prev'],
                                           [-0.1, 0.25, 'prev'],
                                           [0, 0.5, 'prev'],
                                           # [-0.15, 0.3, 'prev']
                                           ]
fo.add_component(fo_bar_2)

# %%%% Fanout bar 3
fo_bar_3 = FanOutLine('barrier_3', 0, collection, fo_points)

fo_bar_3.fo_direction = 'top'
fo_bar_3.n_fanout = 9

fo_bar_3.fo_line_fine.fo_width_start = 40e-3
fo_bar_3.fo_line_fine.points_along_path = [[0.1, 0.1, 'start'],
                                           [0.15, 0.05, 'prev'],
                                           [0.1, 0.25, 'prev'],
                                           [0, 0.3, 'prev'],
                                           [0.15, 0.3, 'prev']
                                           ]
fo.add_component(fo_bar_3)

# %%% Fanout sensor
# %%%% Fanout sensor plunger top
fo_sens_pl_top = FanOutLine('sensor_top_plunger', 0, collection, fo_points)

fo_sens_pl_top.fo_direction = 'top'
fo_sens_pl_top.n_fanout = 5

fo_sens_pl_top.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_top.fo_line_fine.points_along_path = [[0, 0.8, 'start'],
                                                 # [0, pl_ver.diameter/2, 'prev'],
                                                 # [-0.2, 0.5, 'prev']
                                                 ]
fo.add_component(fo_sens_pl_top)

# %%%% Fanout bar sens source
fo_sens_top_sou = FanOutLine(
    'sensor_top_barrier_source', 0, collection, fo_points)

fo_sens_top_sou.fo_direction = 'top'
fo_sens_top_sou.n_fanout = 4

fo_sens_top_sou.fo_line_fine.fo_width_start = 40e-3
fo_sens_top_sou.fo_line_fine.points_along_path = [[0, 0.1, 'start'],
                                                  [0, 0.8, 'prev'],
                                                  # [0, 0.3, 'prev'],
                                                  # [-0.15, 0.3, 'prev']
                                                  ]
fo.add_component(fo_sens_top_sou)

# %%%% Fanout bar sens drain
fo_sens_top_dra = FanOutLine('sensor_top_barrier_drain', 0,
                             collection, fo_points)

fo_sens_top_dra.fo_direction = 'top'
fo_sens_top_dra.n_fanout = 6

fo_sens_top_dra.fo_line_fine.fo_width_start = 40e-3
fo_sens_top_dra.fo_line_fine.points_along_path = [[0, 0.7, 'start'],
                                                  # [0, 0.8, 'prev'],
                                                  # [0, 0.3, 'prev'],
                                                  # [-0.15, 0.3, 'prev']
                                                  ]
fo.add_component(fo_sens_top_dra)


# %%%% Fanout bar sens top sep
fo_sens_top_sep = FanOutLine('sensor_top_barrier_seperation', 0,
                             collection, fo_points)

fo_sens_top_sep.fo_direction = 'top'
fo_sens_top_sep.n_fanout = 2

fo_sens_top_sep.fo_line_fine.fo_width_start = 40e-3
fo_sens_top_sep.fo_line_fine.points_along_path = [[-0.2, 0, 'start'],
                                                  [-0.1, 0.1, 'prev'],
                                                  [0, 0.6, 'prev'],
                                                  # [-0.15, 0.3, 'prev']
                                                  ]
fo.add_component(fo_sens_top_sep)

# %%% Fanout clavier
# %%%% Fanout clavier gate top 1
fo_clav_gate_1 = FanOutLine('clavier_gate_0', 0, collection, fo_points)

fo_clav_gate_1.fo_direction = 'left'
fo_clav_gate_1.n_fanout = 1

fo_clav_gate_1.fo_line_fine.fo_width_start = 100e-3
fo_clav_gate_1.fo_line_fine.points_along_path = [[0, 1, 'start'],
                                                 [-0.4, 0.7, 'prev']]
fo.add_component(fo_clav_gate_1)

# %%%% Fanout clavier gate top 2
fo_clav_gate_1 = FanOutLine('clavier_gate_2', 0, collection, fo_points)

fo_clav_gate_1.fo_direction = 'left'
fo_clav_gate_1.n_fanout = 0

fo_clav_gate_1.fo_line_fine.fo_width_start = 100e-3
fo_clav_gate_1.fo_line_fine.points_along_path = [[0, 2, 'start']]
fo.add_component(fo_clav_gate_1)

# %%%% Fanout clavier gate bottom 1
fo_clav_gate_1 = FanOutLine('clavier_gate_1', 0, collection, fo_points)

fo_clav_gate_1.fo_direction = 'bottom'
fo_clav_gate_1.n_fanout = 0

fo_clav_gate_1.fo_line_fine.fo_width_start = 100e-3
fo_clav_gate_1.fo_line_fine.points_along_path = [[0, -1, 'start'],
                                                 [-0.4, -0.7, 'prev']]
fo.add_component(fo_clav_gate_1)

# %%%% Fanout clavier gate bottom 2
fo_clav_gate_1 = FanOutLine('clavier_gate_3', 0, collection, fo_points)

fo_clav_gate_1.fo_direction = 'bottom'
fo_clav_gate_1.n_fanout = 4

fo_clav_gate_1.fo_line_fine.fo_width_start = 100e-3
fo_clav_gate_1.fo_line_fine.points_along_path = [[0, -2, 'start']]
fo.add_component(fo_clav_gate_1)

# %% Screening fanout

# %%% Plunger 0 screening
screen_pl_0 = ScreeningGate('screening_gate_pl_0', collection)

screen_pl_0.screen('plunger', 0,
                   [0.1, 0.3, 0.5], [60e-3, 60e-3, (80e-3, 30e-3)])
screen_pl_0.layer = screening_layer
screen_pl_0.fo_contact_width = 50e-3
screen_pl_0.fo_contact_direction = 0

screen_pl_0_qda = qda.add_component()
screen_pl_0_qda.component = screen_pl_0

# %%% Plunger 1 screening

screen_pl_1 = ScreeningGate('screening_gate_pl_1', collection)

screen_pl_1.screen('plunger', 1, [0.1, 0.4], [50e-3, 50e-3])
screen_pl_1.layer = screening_layer


screen_pl_1_qda = qda.add_component()
screen_pl_1_qda.component = screen_pl_1

# %%% Plunger 2 screening

screen_pl_2 = ScreeningGate('screening_gate_pl_2', collection)

screen_pl_2.screen('plunger', 2, [0.1, 0.4], [50e-3, 50e-3])
screen_pl_2.layer = screening_layer

screen_pl_2_qda = qda.add_component()
screen_pl_2_qda.component = screen_pl_2

# %% Fanout plunger
# First we have to update the quantum dot array with the screening gates
qda.build()

# %%% Fanout screening plunger sensor top

fo_screen_pl_0_top = FanOutLine('screening_gate_pl_0', 0,
                                collection, fo_points)

fo_screen_pl_0_top.fo_direction = 'top'
fo_screen_pl_0_top.n_fanout = 7

fo_screen_pl_0_top.fo_line_fine.fo_width_start = screen_pl_0.fo_contact_width
fo_screen_pl_0_top.fo_line_fine.points_along_path = [[0.0, 0.1, 'start'],
                                                     [0.05, 0.8, 'start'],
                                                     ]
fo.add_component(fo_screen_pl_0_top)


# %%% Build and Add fanout to qda

fo_qda = qda.add_component()
fo_qda.component = fo

# %% Build and save

qda.build()
layout_path = "..\example_notebooks\layout_files\chip_layout.gds"
layout = qda.add_chip_layout(layout_path)
qda.save_as_gds('..\example_notebooks\example_clavier_device.gds')
