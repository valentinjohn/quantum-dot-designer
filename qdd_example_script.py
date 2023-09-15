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
qda_components = qdd.QuantumDotArrayComponents(qda_elements)
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
pl_ver.asymx = 0.8
pl_ver.diameter = 150e-3
pl_ver.layer = 21

pl_hor = qda_elements.add_plunger('plunger_horizontally_elongated')
pl_hor.asymx = 1.1
pl_hor.diameter = 100e-3
pl_hor.layer = 21

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

# %% define sensor

qda_elements.spacing_sep = 60e-3

sensor_top = qda_components.add_sensor('sensor_top')

# for the positioning of the sensor element you can either indicate 'top',
# 'bottom', 'top-right', etc ...., or you indiacte the angle with 'top'
# corresponding to 0, 'top-right' corresponding to np.pi/4, etc ...
sensor_top.sep_pos = 'bottom'
# sensor_top.sep_pos_angle = 1/8*np.pi
# sensor_top.source_pos = 'right'
sensor_top.source_pos_angle = 3/8*np.pi
# sensor_top.drain_pos = 'left'  # 1/8*np.pi
sensor_top.drain_pos_angle = -3/8*np.pi

# the orientation of the sensor elements is indicated by the direction in
# which the end of the element points to, This can be either clockwise or
# counterclockwise
sensor_top.bar_drain_end = 'counterclockwise'
sensor_top.bar_source_end = 'clockwise'
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
sensor_top.source.contact_length = 70e-3
sensor_top.drain.layer = ohmic_layer
sensor_top.drain.contact_length = 70e-3

sensor_bottom = qda_components.add_copy(sensor_top, 'sensor_bottom')
sensor_bottom.sep_pos = 'top'
sensor_bottom.source_pos = 'left'
sensor_bottom.drain_pos = 'right'

sensor_right = qda_components.add_copy(sensor_top, 'sensor_right')
sensor_right.sep_pos = 'left'
sensor_right.source_pos = 'bottom'
sensor_right.drain_pos = 'top'

sensor_left = qda_components.add_copy(sensor_top, 'sensor_left')
sensor_left.sep_pos = 'right'
sensor_left.source_pos = 'top'
sensor_left.drain_pos = 'bottom'

# sensor_top.plot(build=True)

# %% add sensor to unit cell

unit_cell_ylim = (qda.spacing_qd_diag/2 + pl_ver.diameter/2 * pl_ver.asymy)

sensor_pos_uniy = (unit_cell_ylim+sensor_top.gap_sep +
                   sensor_top.plunger.diameter/2 *
                   sensor_top.plunger.asymy)

sensor_pos_unix = qda.spacing_qd_diag/2

uc_st = unit_cell.add_component()
uc_st.component = sensor_top
uc_st.center = (sensor_pos_unix, sensor_pos_uniy)

uc_sb = unit_cell.add_component()
uc_sb.component = sensor_bottom
uc_sb.center = (-sensor_pos_unix, -sensor_pos_uniy)

# %% Add unit cell to main cell

uc_unitcell = qda.add_component()
uc_unitcell.component = unit_cell
uc_unitcell.center = (0, 0)
uc_unitcell.rows = 1
uc_unitcell.columns = 1
uc_unitcell.spacing = (2*qda.spacing_qd_diag,
                       qda.spacing_qd_diag)

# %% Sensor positions

unit_cell_xlim = (qda.spacing_qd_diag + pl_hor.diameter/2 * pl_hor.asymx)

sensor_pos_x = (unit_cell_xlim +
                sensor_right.gap_sep +
                sensor_right.plunger.diameter/2 *
                sensor_right.plunger.asymx)

uc_sr = qda.add_component()
uc_sr.component = sensor_right
uc_sr.center = (sensor_pos_x, 0)

uc_sl = qda.add_component()
uc_sl.component = sensor_left
uc_sl.center = (-sensor_pos_x, 0)

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


# %%% Fanout plunger ver 0
fo_pl_ver_0 = qda_components.add_fo_line('plunger_vertically_elongated', 0)

fo_pl_ver_0.fo_direction = 'top'
fo_pl_ver_0.n_fanout = 0
fo_pl_ver_0.fo_points = fo_points

fo_pl_ver_0.fo_line_fine.fo_width_start = 40e-3
fo_pl_ver_0.fo_line_fine.points_along_path = [[-0.15, 0.3, 'start'],
                                              [-0.4, 0.7, 'start']]
fo.add_component(fo_pl_ver_0)

# %%% Fanout plunger ver 1
fo_pl_ver_1 = qda_components.add_fo_line('plunger_vertically_elongated', 1)

fo_pl_ver_1.fo_direction = 'top'
fo_pl_ver_1.n_fanout = 2
fo_pl_ver_1.fo_points = fo_points

fo_pl_ver_1.fo_line_fine.fo_width_start = 40e-3
fo_pl_ver_1.fo_line_fine.points_along_path = [[-qda.spacing_qd_diag/2, 0, 'start'],
                                              [0, pl_ver.diameter/2, 'prev'],
                                              [-0.2, 0.5, 'prev']]
fo.add_component(fo_pl_ver_1)

# %%% Fanout bar 45deg 0
fo_bar_45deg_0 = qda_components.add_fo_line('barrier_45deg_rotated', 0)

fo_bar_45deg_0.fo_direction = 'top'
fo_bar_45deg_0.n_fanout = 1
fo_bar_45deg_0.fo_points = fo_points

fo_bar_45deg_0.fo_line_fine.fo_width_start = 40e-3
fo_bar_45deg_0.fo_line_fine.points_along_path = [[0.01, 0.015, 'start'],
                                                 [0.01, 0.05, 'start'],
                                                 [0, 0.08, 'prev'],
                                                 [-0.15, 0.3, 'prev']]
fo.add_component(fo_bar_45deg_0)

# %%% Fanout sensor plunger top
fo_sens_pl_top = qda_components.add_fo_line('sensor_top_plunger', 0)

fo_sens_pl_top.fo_direction = 'top'
fo_sens_pl_top.n_fanout = 4
fo_sens_pl_top.fo_points = fo_points

fo_sens_pl_top.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_top.fo_line_fine.points_along_path = [[0, 0.8, 'start', 40e-3],
                                                 # [0, pl_ver.diameter/2, 'prev'],
                                                 # [-0.2, 0.5, 'prev']
                                                 ]
fo.add_component(fo_sens_pl_top)

# %%% Fanout sensor source top
fo_sens_top_source = qda_components.add_fo_line('sensor_top_source', 0)

fo_sens_top_source.fo_direction = 'top'
fo_sens_top_source.n_fanout = 6
fo_sens_top_source.fo_points = fo_points

fo_sens_top_source.fo_line_fine.points_along_path = [[0.25, 0.8, 'start'],
                                                     # [0, pl_ver.diameter/2, 'prev'],
                                                     # [-0.2, 0.5, 'prev']
                                                     ]
fo.add_component(fo_sens_top_source)

# %%% Fanout sensor drain top
fo_sens_top_drain = qda_components.add_fo_line('sensor_top_drain', 0)

fo_sens_top_drain.fo_direction = 'top'
fo_sens_top_drain.n_fanout = 3
fo_sens_top_drain.fo_points = fo_points

fo_sens_top_drain.fo_line_fine.points_along_path = [[0.25, 0.8, 'start'],
                                                    # [0, pl_ver.diameter/2, 'prev'],
                                                    # [-0.2, 0.5, 'prev']
                                                    ]
fo.add_component(fo_sens_top_drain)

# %%% Fanout sensor plunger bottom
fo_sens_pl_bottom = qda_components.add_fo_line('sensor_bottom_plunger', 0)

fo_sens_pl_bottom.fo_direction = 'bottom'
fo_sens_pl_bottom.n_fanout = 2
fo_sens_pl_bottom.fo_points = fo_points

fo_sens_pl_bottom.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_bottom.fo_line_fine.points_along_path = [[0, -0.8, 'start'],
                                                    # [0, pl_ver.diameter/2, 'prev'],
                                                    # [-0.2, 0.5, 'prev']
                                                    ]
fo.add_component(fo_sens_pl_bottom)

# %%% Fanout sensor plunger right
fo_sens_pl_right = qda_components.add_fo_line('sensor_right_plunger', 0)

fo_sens_pl_right.fo_direction = 'right'
fo_sens_pl_right.n_fanout = 2
fo_sens_pl_right.fo_points = fo_points

fo_sens_pl_right.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_right.fo_line_fine.points_along_path = [[0.8, 0, 'start'],
                                                   [0.9, 0, 'start'],
                                                   # [0, pl_ver.diameter/2, 'prev'],
                                                   # [-0.2, 0.5, 'prev']
                                                   ]
fo.add_component(fo_sens_pl_right)

# %%% Fanout sensor source right

fo_sens_right_source = qda_components.add_fo_line('sensor_right_source', 0)

fo_sens_right_source.fo_direction = 'right'
fo_sens_right_source.n_fanout = 4
fo_sens_right_source.fo_points = fo_points

# fo_sens_right_source.fo_line_fine.points_along_path = [[0.8, 0.25, 'start'],
#                                                        # [0, pl_ver.diameter/2, 'prev'],
#                                                        # [-0.2, 0.5, 'prev']
#                                                        ]
fo.add_component(fo_sens_right_source)

# %%% Fanout sensor drain right

fo_sens_right_drain = qda_components.add_fo_line('sensor_right_drain', 0)

fo_sens_right_drain.fo_direction = 'right'
fo_sens_right_drain.n_fanout = 0
fo_sens_right_drain.fo_points = fo_points

# fo_sens_right_drain.fo_line_fine.points_along_path = [[0.8, 0.25, 'start'],
#                                                        # [0, pl_ver.diameter/2, 'prev'],
#                                                        # [-0.2, 0.5, 'prev']
#                                                        ]
fo.add_component(fo_sens_right_drain)

# %%% Fanout sensor barrier source right

fo_sens_right_barrier_source = qda_components.add_fo_line(
    'sensor_right_barrier_source', 0)

fo_sens_right_barrier_source.fo_direction = 'right'
fo_sens_right_barrier_source.n_fanout = 3
fo_sens_right_barrier_source.fo_points = fo_points

fo_sens_right_barrier_source.fo_line_fine.fo_width_start = 30e-3

fo_sens_right_barrier_source.fo_line_fine.points_along_path = [[0.8, 0, 'start'],
                                                               # [0, pl_ver.diameter/2, 'prev'],
                                                               # [-0.2, 0.5, 'prev']
                                                               ]
fo.add_component(fo_sens_right_barrier_source)

# %%% Fanout sensor drain right
# fo_sens_right_drain = qda_components.add_fo_line('sensor_right_drain', 0)

# fo_sens_right_drain.fo_direction = 'right'
# fo_sens_right_drain.n_fanout = 2
# fo_sens_right_drain.fo_points = fo_points

# fo_sens_right_drain.fo_line_fine.points_along_path = [[0.8, 0.25, 'start'],
#                                                       # [0, pl_ver.diameter/2, 'prev'],
#                                                       # [-0.2, 0.5, 'prev']
#                                                       ]
# fo.add_component(fo_sens_right_drain)

# %%% Fanout sensor plunger left
fo_sens_pl_left = qda_components.add_fo_line('sensor_left_plunger', 0)

fo_sens_pl_left.fo_direction = 'left'
fo_sens_pl_left.n_fanout = 2
fo_sens_pl_left.fo_points = fo_points

fo_sens_pl_left.fo_line_fine.fo_width_start = 40e-3
fo_sens_pl_left.fo_line_fine.points_along_path = [[-0.8, 0, 'start'],
                                                  # [0, pl_ver.diameter/2, 'prev'],
                                                  # [-0.2, 0.5, 'prev']
                                                  ]
fo.add_component(fo_sens_pl_left)

# %% Build and Add fanout to qda

fo_qda = qda.add_component()
fo_qda.component = fo

# %% Build and save

qda.build()
layout = qda.add_chip_layout()
qda.save_as_gds('qdd_test_design.gds')
