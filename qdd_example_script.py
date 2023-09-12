# -*- coding: utf-8 -*-
"""
Created on Tue May  9 08:18:02 2023

@author: vjohn
"""
# %% import

import numpy as np
import QuantumDotDesigner as qdd

# %% Init

qda_elements = qdd.QuantumDotArrayElements()
unit_cell = qdd.UnitCell('unit_cell')
qda = qdd.QuantumDotArray()

# %% Layers

ohmic_layer = 4
barrier_layer = 5
screening_layer = 9
barrier_source_layer = 5
barrier_drain_layer = 5
plunger_layer = 7

clav_gate_1_layer = 5
clav_gate_2_layer = 7

# %% define plungers

pl_ver = qda_elements.add_plunger('plunger_vertically_elongated')
pl_ver.asym = 0.8
pl_ver.diameter = 150e-3
pl_ver.layer = 21

pl_hor = qda_elements.add_plunger('plunger_horizontally_elongated')
pl_hor.asym = 1.1
pl_hor.diameter = 100e-3
pl_hor.layer = 21

pl_ver.build()
pl_hor.build()

# %% define barriers

bar_45deg = qda_elements.add_barrier('barrier_45deg_rotated')
bar_45deg.width = 40e-3
bar_45deg.length = 70e-3
bar_45deg.layer = 5
bar_45deg.rotate = 1/4*np.pi + 1/16*np.pi

bar_135deg = qda_elements.add_copy(bar_45deg, 'barrier_135deg_rotated')
bar_135deg.rotate = 3/4*np.pi - 1/16*np.pi

bar_225deg = qda_elements.add_copy(bar_45deg, 'barrier_225deg_rotated')
bar_225deg.rotate = -3/4*np.pi + 1/16*np.pi

bar_315deg = qda_elements.add_copy(bar_45deg, 'barrier_315deg_rotated')
bar_315deg.rotate = -1/4*np.pi - 1/16*np.pi

bar_45deg.build()
bar_135deg.build()
bar_225deg.build()
bar_315deg.build()

# %% define unit cell

unit_cell.spacing_qd = 200e-3

uc_pl_ver = unit_cell.add_component()
uc_pl_hor = unit_cell.add_component()

uc_bar_45deg = unit_cell.add_component()
uc_bar_135deg = unit_cell.add_component()
uc_bar_225deg = unit_cell.add_component()
uc_bar_315deg = unit_cell.add_component()

uc_pl_ver.component = pl_ver
uc_pl_ver.center = (0, 0)
uc_pl_ver.rows = 2
uc_pl_ver.columns = 2
uc_pl_ver.spacing = (qda.spacing_qd_diag,
                     qda.spacing_qd_diag)

uc_pl_hor.component = pl_hor
uc_pl_hor.center = (0, 0)
uc_pl_hor.columns = 3
uc_pl_hor.spacing = (qda.spacing_qd_diag,
                     qda.spacing_qd_diag)

uc_bar_45deg.component = bar_45deg
uc_bar_45deg.center = (1.15*qda.spacing_qd_diag/4, qda.spacing_qd_diag/4)
uc_bar_45deg.columns = 2
uc_bar_45deg.spacing = (qda.spacing_qd_diag, qda.spacing_qd_diag)

uc_bar_135deg.component = bar_135deg
uc_bar_135deg.center = (-1.15*qda.spacing_qd_diag/4, qda.spacing_qd_diag/4)
uc_bar_135deg.columns = 2
uc_bar_135deg.spacing = (qda.spacing_qd_diag, qda.spacing_qd_diag)

uc_bar_225deg.component = bar_225deg
uc_bar_225deg.center = (-1.15*qda.spacing_qd_diag/4, -qda.spacing_qd_diag/4)
uc_bar_225deg.columns = 2
uc_bar_225deg.spacing = (qda.spacing_qd_diag, qda.spacing_qd_diag)

uc_bar_315deg.component = bar_315deg
uc_bar_315deg.center = (+1.15*qda.spacing_qd_diag/4, -qda.spacing_qd_diag/4)
uc_bar_315deg.columns = 2
uc_bar_315deg.spacing = (qda.spacing_qd_diag, qda.spacing_qd_diag)

uc_pl_ver.build()
uc_pl_hor.build()

uc_bar_135deg.build()
uc_bar_45deg.build()
uc_bar_315deg.build()
uc_bar_225deg.build()

# %% define sensor

qda_elements.spacing_sep = 60e-3

sensor_top = qda_elements.add_sensor('sensor_top')

sensor_top.sep_pos = 5/8*np.pi
sensor_top.source_pos = -7/8*np.pi
sensor_top.drain_pos = 3/8*np.pi
sensor_top.bar_drain_end = 'clockwise'
sensor_top.bar_source_end = 'counter-clockwise'
sensor_top.bar_sep_end = 'clockwise'

sensor_top.gap_sep = 60e-3
sensor_top.gap_ohmic_pl = 40e-3

sensor_top.plunger.diameter = 160e-3
sensor_top.plunger.layer = 21

sensor_top.barrier_source.width = 30e-3
sensor_top.barrier_source.length = 70e-3
sensor_top.barrier_source.layer = 5

sensor_top.barrier_drain.width = 30e-3
sensor_top.barrier_drain.length = 70e-3
sensor_top.barrier_drain.layer = 5

sensor_top.barrier_sep.width = 50e-3
sensor_top.barrier_sep.length = 60e-3
sensor_top.barrier_sep.layer = 5

sensor_top.source.layer = ohmic_layer
sensor_top.drain.layer = ohmic_layer


sensor_bottom = qda_elements.add_copy(sensor_top, 'sensor_bottom')
sensor_bottom.sep_pos = 'top'
sensor_bottom.source.ohmic_pos = 'left'
sensor_bottom.source.sensor_pos = 'bottom'
sensor_bottom.drain.ohmic_pos = 'right'
sensor_bottom.drain.sensor_pos = 'bottom'

sensor_right = qda_elements.add_copy(sensor_top, 'sensor_right')
sensor_right.sep_pos = 'left'
# sensor_right.sensor_pos = 'right'
sensor_right.source.sensor_pos = 'right'  # redundant
sensor_right.drain.sensor_pos = 'right'  # redundant
sensor_right.source_pos = 'top'
sensor_right.drain_pos = 'bottom'
sensor_right.source.ohmic_pos = sensor_right.source_pos  # redundant
sensor_right.drain.ohmic_pos = sensor_right.drain_pos  # redundant

sensor_left = qda_elements.add_copy(sensor_top, 'sensor_left')
sensor_left.sep_pos = 'right'
# sensor_left.sensor_pos = 'right'
sensor_left.source.sensor_pos = 'left'  # redundant
sensor_left.drain.sensor_pos = 'left'  # redundant
sensor_left.source_pos = 'top'
sensor_left.drain_pos = 'bottom'
sensor_left.source.ohmic_pos = sensor_right.source_pos  # redundant
sensor_left.drain.ohmic_pos = sensor_right.drain_pos  # redundant

sensor_top.build()
sensor_bottom.build()
sensor_right.build()
sensor_left.build()

sensor_top.plot()

# %% add sensor to unit cell

unit_cell.ylim = (qda.spacing_qd_diag/2 +
                  pl_ver.diameter/2 *
                  pl_ver._asymy)

sensor_pos_uniy = (unit_cell.ylim+sensor_top.gap_sep +
                   sensor_top.plunger.diameter/2 *
                   sensor_top.plunger._asymy)

sensor_pos_unix = qda.spacing_qd_diag/2

uc_st = unit_cell.add_component()
uc_st.component = sensor_top
uc_st.center = (sensor_pos_unix, sensor_pos_uniy)
uc_st.build()

uc_sb = unit_cell.add_component()
uc_sb.component = sensor_bottom
uc_sb.center = (-sensor_pos_unix, -sensor_pos_uniy)
uc_sb.build()

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

# %% Sensor positions

sensor_pos_x = (uc_unitcell.xlim[1] +
                sensor_right.gap_sep +
                sensor_right.plunger.diameter/2 *
                sensor_right.plunger._asymx)

uc_sr = qda.add_component()
uc_sr.component = sensor_right
uc_sr.center = (sensor_pos_x, 0)
uc_sr.build()

uc_sl = qda.add_component()
uc_sl.component = sensor_left
uc_sl.center = (-sensor_pos_x, 0)
uc_sl.build()

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


# %% Fanout plunger ver 0
fo_pl_ver_0 = qda_elements.add_fo_line('plunger_vertically_elongated', 0)

fo_pl_ver_0.fo_direction = 'top'
fo_pl_ver_0.n_fanout = 0
fo_pl_ver_0.fo = fo

fo_pl_ver_0.fo_line_fine.fo_width_start = 40e-3
fo_pl_ver_0.fo_line_fine.points_along_path = [[-0.15, 0.3, 'start'],
                                              [-0.4, 0.7, 'start']]
fo_pl_ver_0.build()
fo.add_component(fo_pl_ver_0)

# %% Fanout plunger ver 1
fo_pl_ver_1 = qda_elements.add_fo_line('plunger_vertically_elongated', 1)

fo_pl_ver_1.fo_direction = 'top'
fo_pl_ver_1.n_fanout = 2
fo_pl_ver_1.fo = fo

fo_pl_ver_1.fo_line_fine.fo_width_start = 40e-3
fo_pl_ver_1.fo_line_fine.points_along_path = [[-qda.spacing_qd_diag/2, 0, 'start'],
                                              [0, pl_ver.diameter/2, 'prev'],
                                              [-0.2, 0.5, 'prev']]
fo_pl_ver_1.build()
fo.add_component(fo_pl_ver_1)

# %% Fanout bar 45deg 0
fo_bar_45deg_0 = qda_elements.add_fo_line('barrier_45deg_rotated', 0)

fo_bar_45deg_0.fo_direction = 'top'
fo_bar_45deg_0.n_fanout = 1
fo_bar_45deg_0.fo = fo

fo_bar_45deg_0.fo_line_fine.fo_width_start = 40e-3
fo_bar_45deg_0.fo_line_fine.points_along_path = [[0.01, 0.015, 'start'],
                                                 [0.01, 0.05, 'start'],
                                                 [0, 0.08, 'prev'],
                                                 [-0.15, 0.3, 'prev']]
fo_bar_45deg_0.build()
fo.add_component(fo_bar_45deg_0)

# %% Fanout sensor plunger top
fo_sens_pl_top = qda_elements.add_fo_line('sensor_top_plunger', 0)

fo_sens_pl_top.fo_direction = 'top'
fo_sens_pl_top.n_fanout = 4
fo_sens_pl_top.fo = fo

fo_sens_pl_top.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_top.fo_line_fine.points_along_path = [[0, 0.8, 'start'],
                                                 # [0, pl_ver.diameter/2, 'prev'],
                                                 # [-0.2, 0.5, 'prev']
                                                 ]
fo_sens_pl_top.build()
fo.add_component(fo_sens_pl_top)

# %% Fanout sensor plunger bottom
fo_sens_pl_bottom = qda_elements.add_fo_line('sensor_bottom_plunger', 0)

fo_sens_pl_bottom.fo_direction = 'bottom'
fo_sens_pl_bottom.n_fanout = 2
fo_sens_pl_bottom.fo = fo

fo_sens_pl_bottom.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_bottom.fo_line_fine.points_along_path = [[0, -0.8, 'start'],
                                                    # [0, pl_ver.diameter/2, 'prev'],
                                                    # [-0.2, 0.5, 'prev']
                                                    ]
fo_sens_pl_bottom.build()
fo.add_component(fo_sens_pl_bottom)

# %% Fanout sensor plunger right
fo_sens_pl_right = qda_elements.add_fo_line('sensor_right_plunger', 0)

fo_sens_pl_right.fo_direction = 'right'
fo_sens_pl_right.n_fanout = 2
fo_sens_pl_right.fo = fo

fo_sens_pl_right.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_right.fo_line_fine.points_along_path = [[0.8, 0, 'start'],
                                                   [0.9, 0, 'start'],
                                                   # [0, pl_ver.diameter/2, 'prev'],
                                                   # [-0.2, 0.5, 'prev']
                                                   ]
fo_sens_pl_right.build()
fo.add_component(fo_sens_pl_right)

# %% Fanout sensor plunger left
fo_sens_pl_left = qda_elements.add_fo_line('sensor_left_plunger', 0)

fo_sens_pl_left.fo_direction = 'left'
fo_sens_pl_left.n_fanout = 2
fo_sens_pl_left.fo = fo

fo_sens_pl_left.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_left.fo_line_fine.points_along_path = [[-0.8, 0, 'start'],
                                                  # [0, pl_ver.diameter/2, 'prev'],
                                                  # [-0.2, 0.5, 'prev']
                                                  ]
fo_sens_pl_left.build()
fo.add_component(fo_sens_pl_left)

# %% Build and Add fanout to qda

fo.build()
fo_qda = qda.add_component()
fo_qda.component = fo
fo_qda.build()

# %% Build and save

qda.build()
layout = qda.add_chip_layout()
qda.save_as_gds('qdd_test_design.gds')
