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

# %% define plungers

pl_ver = qda_elements.add_plunger('plunger_vertically_elongated')
pl_ver.asym = 0.8
pl_ver.diameter = 150
pl_ver.layer = 21

pl_hor = qda_elements.add_plunger('plunger_horizontally_elongated')
pl_hor.asym = 1.1
pl_hor.diameter = 100
pl_hor.layer = 21

pl_ver.build()
pl_hor.build()

# %% define barriers

bar_one_fourths_pi = qda_elements.add_barrier('barrier_one_fourths_pi_rotated')
bar_one_fourths_pi.width = 30
bar_one_fourths_pi.length = 70
bar_one_fourths_pi.layer = 5
bar_one_fourths_pi.rotate = 1/4*np.pi

bar_three_fourths_pi = qda_elements.add_copy(bar_one_fourths_pi,
                                             'barrier_three_fourths_pi_rotated')
bar_three_fourths_pi.rotate = 3/4*np.pi

bar_five_fourths_pi = qda_elements.add_copy(bar_one_fourths_pi,
                                            'barrier_five_fourths_pi_rotated')
bar_five_fourths_pi.rotate = -3/4*np.pi

bar_seven_fourths_pi = qda_elements.add_copy(bar_one_fourths_pi,
                                             'barrier_seven_fourths_pi_rotated')
bar_seven_fourths_pi.rotate = -1/4*np.pi

bar_one_fourths_pi.build()
bar_three_fourths_pi.build()
bar_five_fourths_pi.build()
bar_seven_fourths_pi.build()

# %% define unit cell

unit_cell.spacing_qd = 200

sl_pl_ver = unit_cell.add_component()
sl_pl_hor = unit_cell.add_component()

sl_bar_three_fourths_pi = unit_cell.add_component()
sl_bar_one_fourths_pi = unit_cell.add_component()
sl_bar_seven_fourths_pi = unit_cell.add_component()
sl_bar_five_fourths_pi = unit_cell.add_component()

sl_pl_ver.component = pl_ver
sl_pl_ver.center = (0, 0)
sl_pl_ver.rows = 2
sl_pl_ver.columns = 2
sl_pl_ver.spacing = (qda.spacing_qd_diag,
                     qda.spacing_qd_diag)

sl_pl_hor.component = pl_hor
sl_pl_hor.center = (0, 0)
sl_pl_hor.columns = 3
sl_pl_hor.spacing = (qda.spacing_qd_diag,
                     qda.spacing_qd_diag)

sl_bar_one_fourths_pi.component = bar_one_fourths_pi
sl_bar_one_fourths_pi.center = (+qda.spacing_qd_diag/4,
                                qda.spacing_qd_diag/4)
sl_bar_one_fourths_pi.columns = 2
sl_bar_one_fourths_pi.spacing = (qda.spacing_qd_diag,
                                 qda.spacing_qd_diag)

sl_bar_three_fourths_pi.component = bar_three_fourths_pi
sl_bar_three_fourths_pi.center = (-qda.spacing_qd_diag/4,
                                  qda.spacing_qd_diag/4)
sl_bar_three_fourths_pi.columns = 2
sl_bar_three_fourths_pi.spacing = (qda.spacing_qd_diag,
                                   qda.spacing_qd_diag)

sl_bar_five_fourths_pi.component = bar_five_fourths_pi
sl_bar_five_fourths_pi.center = (-qda.spacing_qd_diag/4,
                                 -qda.spacing_qd_diag/4)
sl_bar_five_fourths_pi.columns = 2
sl_bar_five_fourths_pi.spacing = (qda.spacing_qd_diag,
                                  qda.spacing_qd_diag)

sl_bar_seven_fourths_pi.component = bar_seven_fourths_pi
sl_bar_seven_fourths_pi.center = (+qda.spacing_qd_diag/4,
                                  -qda.spacing_qd_diag/4)
sl_bar_seven_fourths_pi.columns = 2
sl_bar_seven_fourths_pi.spacing = (qda.spacing_qd_diag,
                                   qda.spacing_qd_diag)

sl_pl_ver.build()
sl_pl_hor.build()

sl_bar_three_fourths_pi.build()
sl_bar_one_fourths_pi.build()
sl_bar_seven_fourths_pi.build()
sl_bar_five_fourths_pi.build()

# %% define sensor

# qda_elements.spacing_sep = 60

# sensor_top = qda_elements.add_sensor('sensor_top')
# sensor_top.source_pos = 'left'
# sensor_top.drain_pos = 'right'
# sensor_top.sep_pos = 'bottom'
# sensor_top.gap_sep = 60
# sensor_top.gap_ohmic_pl = 40

# sensor_top.plunger.diameter = 160
# sensor_top.plunger.layer = 21

# sensor_top.barrier_source.width = 30
# sensor_top.barrier_source.length = 70
# sensor_top.barrier_source.layer = 5

# sensor_top.barrier_drain.width = 30
# sensor_top.barrier_drain.length = 70
# sensor_top.barrier_drain.layer = 5

# sensor_top.barrier_sep.width = 50
# sensor_top.barrier_sep.length = 60
# sensor_top.barrier_sep.layer = 5

# sensor_bottom = qda_elements.add_copy(sensor_top, 'sensor_bottom')
# sensor_bottom.source_pos = 'right'
# sensor_bottom.drain_pos = 'left'
# sensor_bottom.sep_pos = 'top'

# sensor_right = qda_elements.add_copy(sensor_top, 'sensor_right')
# sensor_right.source_pos = 'top'
# sensor_right.drain_pos = 'bottom'
# sensor_right.sep_pos = 'left'

# sensor_left = qda_elements.add_copy(sensor_top, 'sensor_left')
# sensor_left.source_pos = 'bottom'
# sensor_left.drain_pos = 'top'
# sensor_left.sep_pos = 'right'

# sensor_top.build()
# sensor_bottom.build()
# sensor_right.build()
# sensor_left.build()

# %% add sensor to unit cell

# unit_cell.ylim = (qda.spacing_qd_diag/2 +
#                   pl_ver.diameter/2 *
#                   pl_ver.get_asym()[1])

# sensor_pos_uniy = (unit_cell.ylim+sensor_top.gap_sep +
#                    sensor_top.plunger.diameter/2 *
#                    sensor_top.plunger.get_asym()[1])

# sensor_pos_unix = qda.spacing_qd_diag/2

# sl_st = unit_cell.add_component('sublattice_sensor_top')
# sl_st.component = sensor_top
# sl_st.center = (sensor_pos_unix, sensor_pos_uniy)
# sl_st.build()

# sl_sb = unit_cell.add_component('sublattice_sensor_bottom')
# sl_sb.component = sensor_bottom
# sl_sb.center = (-sensor_pos_unix, -sensor_pos_uniy)
# sl_sb.build()

# %% Add unit cell to main cell

unit_cell.build()

sl_unitcell = qda.add_component('sublattice_unit_cell')
sl_unitcell.component = unit_cell
sl_unitcell.center = (0, 0)
sl_unitcell.rows = 1
sl_unitcell.columns = 10
sl_unitcell.spacing = (2*qda.spacing_qd_diag,
                       qda.spacing_qd_diag)
sl_unitcell.build()

# %% Sensor positions

# sensor_pos_x = (sl_unitcell.xlim[1] +
#                 sensor_right.gap_sep +
#                 sensor_right.plunger.diameter/2 *
#                 sensor_right.plunger.get_asym()[0])

# sl_sr = qda.add_component('sublattice_sensor_right')
# sl_sr.component = sensor_right
# sl_sr.center = (sensor_pos_x, 0)
# sl_sr.build()

# sl_sl = qda.add_component('sublattice_sensor_left')
# sl_sl.component = sensor_left
# sl_sl.center = (-sensor_pos_x, 0)
# sl_sl.build()

# %% Build and save

qda.build()
qda.save_as_gds('qdd_test_design.gds')


# %%
fog = FanoutGenerator()

(P1, P2, P3, P8, P9, P10) = fog.gates(sl_pl_ver)
(P4, P5, P6, P7) = fog.gates(sl_pl_hor)

P1.points
P1.position
P1.layer
P1.add_fanout()

# #%% Generate plungers and barriers

# P1 = sl_pl_ver.gates[0]
# P2 = sl_pl_ver.gates[1]
# P3 = sl_pl_ver.gates[2]
# P4 = sl_pl_hor.gates[0]
# P5 = sl_pl_hor.gates[1]
# P6 = sl_pl_hor.gates[2]
# P7 = sl_pl_hor.gates[3]
# P8 = sl_pl_ver.gates[3]
# P9 = sl_pl_ver.gates[4]
# P10 = sl_pl_ver.gates[5]


# B1 = sl_bar_three_fourths_pi[0]
# B2 = sl_bar_one_fourths_pi[0]
# B3 = sl_bar_three_fourths_pi[1]
# B4 = sl_bar_one_fourths_pi[1]
# B5 = sl_bar_three_fourths_pi[2]
# B6 = sl_bar_one_fourths_pi[2]
# B7 = sl_bar_seven_fourths_pi[0]
# B8 = sl_bar_five_fourths_pi[0]
# B9 = sl_bar_seven_fourths_pi[1]
# B10 = sl_bar_five_fourths_pi[1]
# B11 = sl_bar_seven_fourths_pi[2]
B12 = sl_bar_five_fourths_pi[2]

# %% define screening gates

TSC = qda_elements.add_screening_gate()
RSC = qda_elements.add_screening_gate()
BSC = qda_elements.add_screening_gate()
LSC = qda_elements.add_screening_gate()

TLSC = qda_elements.add_screening_gate()
TRSC = qda_elements.add_screening_gate()
BRSC = qda_elements.add_screening_gate()
BLSC = qda_elements.add_screening_gate()

TSC.screened_gates = sensor_top
RSC.screened_gates = sensor_right
BSC.screened_gates = sensor_bottom
LSC.screened_gates = sensor_left

TLSC.screened_gates = [P1, P2]
TRSC.screened_gates = [P3, P6, P7]
BRSC.screened_gates = [P9, P10]
BLSC.screened_gates = [P4, P5, P8]

# %% define fine fan-out
fog = FanoutGenerator()

fog.fine_fo_window = (4.5e3, 3e3)
fog.broad_fo_window = (30e3, 30e3)

fog.ohmics = [TOL, TOR, ROT, ROB, BOR, BOL, LOB, LOT]
fog.assign_top_fo = [B2, P2, B3, P3, TOL,
                     TBL, S_T, TSC, TBR, TOR, TS, B4, P6, B5]
fog.assign_right_fo = list()
fog.assign_bottom_fo = list()
fog.assign_left_fo = list()

# %% add fine fan-out

P1_via = P1.add_via()

P1_fo = P1.add_fanout()
P2_fo = P2.add_fanout()
P3_fo = P3.add_fanout()
P4_fo = P4.add_fanout()

P1_via.scale = 0.8
P1_fo.add_relative_points([[0.1, 0.1],
                           [0.2, 0.5],
                           [0.4, 0.8]])


P2_fo.add_relative_points([[0.1, 0.1],
                           [0.2, 0.5],
                           [0.4, 0.8]])

# %% Ideas

# P1 and all other elements should be defined as an object with the following attributes
P1.position
P1.points
P1.raw_def
P1.add_via
P1.add_fanout
