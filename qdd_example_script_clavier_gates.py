# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 16:28:10 2023

@author: vjohn
"""

# %% import

import numpy as np
import QuantumDotDesigner as qdd

# %% Init

qda_elements = qdd.QuantumDotArrayElements()
unit_cell = qdd.UnitCell('unit_cell')
qda = qdd.QuantumDotArray()

# %% define plungers

pl = qda_elements.add_plunger('plunger')
pl.asym = 1
pl.diameter = 130e-3
pl.layer = 3
pl.build()

# %% define barriers

bars_rot = np.pi*(-1/2+np.array([0, -2/3, 2/3]))

bar_1 = qda_elements.add_barrier('barrier_1')
bar_1.width = 40e-3
bar_1.length = 30e-3
bar_1.layer = 1
bar_1.rotate = bars_rot[0]

bar_2 = qda_elements.add_copy(bar_1, 'barrier_2')
bar_2.rotate = bars_rot[1]

bar_3 = qda_elements.add_copy(bar_1, 'barrier_3')
bar_3.rotate = bars_rot[2]

bar_1.build()
bar_2.build()
bar_3.build()

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

uc_pl_tribase.build()
uc_pl_tritop.build()
uc_bar_2.build()
uc_bar_1.build()
uc_bar_3.build()

# %% define sensor

sensor_top = qda_elements.add_sensor('sensor_top')
sensor_top.source_pos = 'left'
sensor_top.drain_pos = 'right'
sensor_top.sep_pos = 'bottom'
sensor_top.gap_sep = 60e-3
sensor_top.gap_ohmic_pl = 40e-3

sensor_top.plunger.diameter = 160e-3
sensor_top.plunger.layer = 21

sensor_top.barrier_source.width = 40e-3
sensor_top.barrier_source.length = 70e-3
sensor_top.barrier_source.layer = 5

sensor_top.barrier_drain.width = 40e-3
sensor_top.barrier_drain.length = 70e-3
sensor_top.barrier_drain.layer = 5

sensor_top.barrier_sep.width = 50e-3
sensor_top.barrier_sep.length = 60e-3
sensor_top.barrier_sep.layer = 5

sensor_top.build()

# %% add sensor to unit cell

unit_cell.ylim = (tritop_height +
                  pl.diameter/2 *
                  pl._asymy)

sensor_pos_uniy = (unit_cell.ylim+sensor_top.gap_sep +
                   sensor_top.plunger.diameter/2 *
                   sensor_top.plunger._asymy)

sensor_pos_unix = 0

uc_st = unit_cell.add_component()
uc_st.component = sensor_top
uc_st.center = (sensor_pos_unix, sensor_pos_uniy)
uc_st.build()

# %% Add unit cell to main cell

unit_cell.build()

uc_unitcell = qda.add_component()
uc_unitcell.component = unit_cell
uc_unitcell.center = (0, 0)
uc_unitcell.rows = 1
uc_unitcell.columns = 1
uc_unitcell.spacing = (2*qda.spacing_qd_diag,
                       qda.spacing_qd_diag)
uc_unitcell.build()

# %% Add clavier gates

clavier = qda_elements.add_clavier('clavier')

clavier.clav_dot_size = 100e-3
clavier.clav_gate_gap = 20e-3
clavier.clav_gap = [.450, .650]
clavier.clav_width = 200e-3
clavier.n_clav_rep = 8

clavier.screen_width = 100e-3
clavier.screen_gap = 0
clavier.screen_layer = 3
clavier.screen_position = 200e-3

clavier.clav_layers = [25, 26]
clavier.x = -(spacing_qd + clavier.screen_length/2)
clavier.y = 0
clavier.fillet = 0.02
clavier.fillet_tolerance = 1e-4

clavier.mirror = False

clavier.build()

clavier_mirrored = qda_elements.add_copy(clavier, 'clavier_mirrored')
clavier_mirrored.mirror = True
clavier_mirrored.x = (spacing_qd + clavier.screen_length/2)
clavier_mirrored.build()

sl_clavier = qda.add_component()
sl_clavier.component = clavier
sl_clavier.build()

sl_clavier_mirrored = qda.add_component()
sl_clavier_mirrored.component = clavier_mirrored
sl_clavier_mirrored.build()

# %% Fanout Generator

qda.build()
fo = qdd.FanoutGenerator('fanout', qda)

fo.fanout_counts = {'top': 10, 'bottom': 10, 'left': 5, 'right': 5}

fo.fo_stages = [(16, 16), (500, 530), (1200, 1200)]
fo.bondpad_position = {'top':  1500, 'bottom': 1500,
                       'left': 1500, 'right': 1500}
fo.fo_widths = [1, 6, 25]
fo.spacings = [2, 80, 250]

fo.bondpad_size = {'top': (110, 400), 'bottom': (110, 400),
                   'left': (400, 110), 'right': (400, 110)}
fo.create_fo_polygons_coarse()

# %% Fanout plunger 0
fo_pl_0 = qda_elements.add_fo_line('plunger', 0,
                                   fo_direction='top', n_fanout=8, fo=fo)

# Fine fan-out
fof_pl_0 = fo_pl_0[0]
fof_pl_0.fo_width_start = 40e-3
fof_pl_0.points_along_path = [[0.25, 0.15, 'start'],
                              [0.15, 0.8, 'prev']]
fof_pl_0.build()
fo.add_component(fof_pl_0)

# Coarse fan-out
foc_pl_0 = fo_pl_0[1]
foc_pl_0.build()
fo.add_component(foc_pl_0)

# %% Fanout plunger 1
fo_pl_1 = qda_elements.add_fo_line('plunger', 1,
                                   fo_direction='bottom', n_fanout=4, fo=fo)

# Fine fan-out
fof_pl_1 = fo_pl_1[0]
fof_pl_1.fo_width_start = 40e-3
fof_pl_1.points_along_path = [[0, -0.1, 'start'],
                              [0, -0.7, 'prev']]
fof_pl_1.build()
fo.add_component(fof_pl_1)

# Coarse fan-out
foc_pl_1 = fo_pl_1[1]
foc_pl_1.build()
fo.add_component(foc_pl_1)

# %% Fanout plunger 2
fo_pl_2 = qda_elements.add_fo_line('plunger', 2,
                                   fo_direction='bottom', n_fanout=6, fo=fo)

# Fine fan-out
fof_pl_2 = fo_pl_2[0]
fof_pl_2.fo_width_start = 40e-3
fof_pl_2.points_along_path = [[0, -0.1, 'start'],
                              [0, -0.7, 'prev']]
fof_pl_2.build()
fo.add_component(fof_pl_2)

# Coarse fan-out
foc_pl_2 = fo_pl_2[1]
foc_pl_2.build()
fo.add_component(foc_pl_2)

# %% Fanout bar 1
fo_bar_1 = qda_elements.add_fo_line('barrier_1', 0,
                                    fo_direction='bottom', n_fanout=5, fo=fo)

# Fine fan-out
fof_bar_1 = fo_bar_1[0]
fof_bar_1.fo_width_start = 40e-3
fof_bar_1.points_along_path = [[0, -0.5, 'start'],
                               # [-0.1, 0.1, 'prev'],
                               # [0, 0.3, 'prev'],
                               # [-0.15, 0.3, 'prev']
                               ]
fof_bar_1.build()
fo.add_component(fof_bar_1)

# Coarse fan-out
foc_bar_1 = fo_bar_1[1]
foc_bar_1.build()
fo.add_component(foc_bar_1)

# %% Fanout bar 2
fo_bar_2 = qda_elements.add_fo_line('barrier_2', 0,
                                    fo_direction='top', n_fanout=1, fo=fo)

# Fine fan-out
fof_bar_2 = fo_bar_2[0]
fof_bar_2.fo_width_start = 40e-3
fof_bar_2.points_along_path = [[-0.1, 0.1, 'start'],
                               [-0.15, 0.05, 'prev'],
                               [-0.1, 0.25, 'prev'],
                               [0, 0.5, 'prev'],
                               # [-0.15, 0.3, 'prev']
                               ]
fof_bar_2.build()
fo.add_component(fof_bar_2)

# Coarse fan-out
foc_bar_2 = fo_bar_2[1]
foc_bar_2.build()
fo.add_component(foc_bar_2)

# %% Fanout bar 3
fo_bar_3 = qda_elements.add_fo_line('barrier_3', 0,
                                    fo_direction='top', n_fanout=9, fo=fo)

# Fine fan-out
fof_bar_3 = fo_bar_3[0]
fof_bar_3.fo_width_start = 40e-3
fof_bar_3.points_along_path = [[0.1, 0.1, 'start'],
                               [0.15, 0.05, 'prev'],
                               [0.1, 0.25, 'prev'],
                               [0, 0.3, 'prev'],
                               [0.15, 0.3, 'prev']
                               ]
fof_bar_3.build()
fo.add_component(fof_bar_3)

# Coarse fan-out
foc_bar_3 = fo_bar_3[1]
foc_bar_3.build()
fo.add_component(foc_bar_3)

# %% Fanout sensor plunger top
fo_sens_pl_top = qda_elements.add_fo_line('sensor_top_plunger', 0,
                                          fo_direction='top', n_fanout=5, fo=fo)

# Fine fan-out
fof_sens_pl_top = fo_sens_pl_top[0]
fof_sens_pl_top.fo_width_start = 40e-3
fof_sens_pl_top.points_along_path = [[0, 0.8, 'start'],
                                     # [0, pl_ver.diameter/2, 'prev'],
                                     # [-0.2, 0.5, 'prev']
                                     ]
fof_sens_pl_top.build()
fo.add_component(fof_sens_pl_top)

# Coarse fan-out
foc_sens_pl_top = fo_sens_pl_top[1]
foc_sens_pl_top.build()
fo.add_component(foc_sens_pl_top)

# %% Fanout bar sens source
fo_sens_top_sou = qda_elements.add_fo_line('sensor_top_barrier_source', 0,
                                           fo_direction='top', n_fanout=4, fo=fo)

# Fine fan-out
fof_sens_top_sou = fo_sens_top_sou[0]
fof_sens_top_sou.fo_width_start = 40e-3
fof_sens_top_sou.points_along_path = [[0, 0.1, 'start'],
                                      [0, 0.8, 'prev'],
                                      # [0, 0.3, 'prev'],
                                      # [-0.15, 0.3, 'prev']
                                      ]
fof_sens_top_sou.build()
fo.add_component(fof_sens_top_sou)

# Coarse fan-out
foc_sens_top_sou = fo_sens_top_sou[1]
foc_sens_top_sou.build()
fo.add_component(foc_sens_top_sou)

# %% Fanout bar sens drain
fo_sens_top_dra = qda_elements.add_fo_line('sensor_top_barrier_drain', 0,
                                           fo_direction='top', n_fanout=6, fo=fo)

# Fine fan-out
fof_sens_top_dra = fo_sens_top_dra[0]
fof_sens_top_dra.fo_width_start = 40e-3
fof_sens_top_dra.points_along_path = [[0, 0.7, 'start'],
                                      # [0, 0.8, 'prev'],
                                      # [0, 0.3, 'prev'],
                                      # [-0.15, 0.3, 'prev']
                                      ]
fof_sens_top_dra.build()
fo.add_component(fof_sens_top_dra)

# Coarse fan-out
foc_sens_top_dra = fo_sens_top_dra[1]
foc_sens_top_dra.build()
fo.add_component(foc_sens_top_dra)

# %% Fanout bar sens top sep
fo_sens_top_sep = qda_elements.add_fo_line('sensor_top_barrier_seperation', 0,
                                           fo_direction='top', n_fanout=2, fo=fo)

# Fine fan-out
fof_sens_top_sep = fo_sens_top_sep[0]
fof_sens_top_sep.fo_width_start = 40e-3
fof_sens_top_sep.points_along_path = [[-0.2, 0, 'start'],
                                      [-0.1, 0.1, 'prev'],
                                      [0, 0.6, 'prev'],
                                      # [-0.15, 0.3, 'prev']
                                      ]
fof_sens_top_sep.calculate_fine_fo()
fof_sens_top_sep.build()
fo.add_component(fof_sens_top_sep)

# Coarse fan-out
foc_sens_top_sep = fo_sens_top_sep[1]
foc_sens_top_sep.build()
fo.add_component(foc_sens_top_sep)

# %% Fanout clavier gate top 1
fo_clav_gate_1 = qda_elements.add_fo_line('clavier_gate_0', 0,
                                          fo_direction='top', n_fanout=2, fo=fo)

# Fine fan-out
fof_clav_gate_1 = fo_clav_gate_1[0]
fof_clav_gate_1.fo_width_start = 100e-3
fof_clav_gate_1.points_along_path = [[0, 1, 'start'],
                                     [-0.4, 0.7, 'prev']]
fof_clav_gate_1.build()
fo.add_component(fof_clav_gate_1)

# Coarse fan-out
foc_clav_gate_1 = fo_clav_gate_1[1]
foc_clav_gate_1.build()
fo.add_component(foc_clav_gate_1)

# %% Fanout clavier gate top 2
fo_clav_gate_1 = qda_elements.add_fo_line('clavier_gate_2', 0,
                                          fo_direction='top', n_fanout=3, fo=fo)

# Fine fan-out
fof_clav_gate_1 = fo_clav_gate_1[0]
fof_clav_gate_1.fo_width_start = 100e-3
fof_clav_gate_1.points_along_path = [[0, 2, 'start']]
fof_clav_gate_1.build()
fo.add_component(fof_clav_gate_1)

# Coarse fan-out
foc_clav_gate_1 = fo_clav_gate_1[1]
foc_clav_gate_1.build()
fo.add_component(foc_clav_gate_1)

# %% Fanout clavier gate bottom 1
fo_clav_gate_1 = qda_elements.add_fo_line('clavier_gate_1', 0,
                                          fo_direction='bottom', n_fanout=0, fo=fo)

# Fine fan-out
fof_clav_gate_1 = fo_clav_gate_1[0]
fof_clav_gate_1.fo_width_start = 100e-3
fof_clav_gate_1.points_along_path = [[0, -1, 'start'],
                                     [-0.4, -0.7, 'prev']]
fof_clav_gate_1.build()
fo.add_component(fof_clav_gate_1)

# Coarse fan-out
foc_clav_gate_1 = fo_clav_gate_1[1]
foc_clav_gate_1.build()
fo.add_component(foc_clav_gate_1)

# %% Fanout clavier gate bottom 2
fo_clav_gate_1 = qda_elements.add_fo_line('clavier_gate_3', 0,
                                          fo_direction='bottom', n_fanout=4, fo=fo)

# Fine fan-out
fof_clav_gate_1 = fo_clav_gate_1[0]
fof_clav_gate_1.fo_width_start = 100e-3
fof_clav_gate_1.points_along_path = [[0, -2, 'start']]
fof_clav_gate_1.build()
fo.add_component(fof_clav_gate_1)

# Coarse fan-out
foc_clav_gate_1 = fo_clav_gate_1[1]
foc_clav_gate_1.build()
fo.add_component(foc_clav_gate_1)

# %% Build and Add fanout to qda

fo.build()
fo_qda = qda.add_component()
fo_qda.component = fo
fo_qda.build()

# %% Build and save

qda.build()
layout = qda.add_chip_layout()
qda.save_as_gds('clavier_gate_design.gds')
