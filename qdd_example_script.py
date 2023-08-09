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
unit_cell = qdd.UnitCell()
qda = qdd.QuantumDotArray()

# %% define plungers

pl_1 = qda_elements.add_plunger('plunger_1')
pl_2 = qda_elements.add_plunger('plunger_2')

pl_1.diameter = 150
pl_1.asym = 0.8
pl_1.layer = 21

pl_2.diameter = 100
pl_2.asym = 1.1
pl_2.layer = 21

pl_1.build()
pl_2.build()

# %% define barriers

bar_qda_1 = qda_elements.add_barrier('barrier_array_1')
bar_qda_1.width = 30
bar_qda_1.length = 70
bar_qda_1.layer = 5
bar_qda_1.rotate = 1/4*np.pi

bar_qda_2 = qda_elements.add_copy(bar_qda_1, 'barrier_array_2')
bar_qda_2.rotate = 3/4*np.pi

bar_qda_3 = qda_elements.add_copy(bar_qda_1, 'barrier_array_3')
bar_qda_3.rotate = -1/4*np.pi

bar_qda_4 = qda_elements.add_copy(bar_qda_1, 'barrier_array_4')
bar_qda_4.rotate = -3/4*np.pi

bar_qda_1.build()
bar_qda_2.build()
bar_qda_3.build()
bar_qda_4.build()

# %% define unit cell

unit_cell.spacing_qd = 200

sl_pl1 = unit_cell.add_component('sublattice_plunger_1')
sl_pl2 = unit_cell.add_component('sublattice_plunger_2')

sl_bar1 = unit_cell.add_component('sublattice_barrier_1')
sl_bar2 = unit_cell.add_component('sublattice_barrier_2')
sl_bar3 = unit_cell.add_component('sublattice_barrier_3')
sl_bar4 = unit_cell.add_component('sublattice_barrier_4')

sl_pl1.component = pl_1
sl_pl1.center = (0, 0)
sl_pl1.rows = 2
sl_pl1.columns = 2
sl_pl1.spacing = (qda.spacing_qd_diag,
                  qda.spacing_qd_diag)

sl_pl2.component = pl_2
sl_pl2.center = (0, 0)
sl_pl2.columns = 3
sl_pl2.spacing = (qda.spacing_qd_diag,
                  qda.spacing_qd_diag)

sl_bar1.component = bar_qda_2
sl_bar1.center = (-qda.spacing_qd_diag/4, qda.spacing_qd_diag/4)
sl_bar1.columns = 2
sl_bar1.spacing = (qda.spacing_qd_diag,
                   qda.spacing_qd_diag)

sl_bar2.component = bar_qda_1
sl_bar2.center = (+qda.spacing_qd_diag/4, qda.spacing_qd_diag/4)
sl_bar2.columns = 2
sl_bar2.spacing = (qda.spacing_qd_diag,
                   qda.spacing_qd_diag)

sl_bar3.component = bar_qda_3
sl_bar3.center = (+qda.spacing_qd_diag/4, -qda.spacing_qd_diag/4)
sl_bar3.columns = 2
sl_bar3.spacing = (qda.spacing_qd_diag,
                   qda.spacing_qd_diag)

sl_bar4.component = bar_qda_4
sl_bar4.center = (-qda.spacing_qd_diag/4, -qda.spacing_qd_diag/4)
sl_bar4.columns = 2
sl_bar4.spacing = (qda.spacing_qd_diag,
                   qda.spacing_qd_diag)

sl_pl1.build()
sl_pl2.build()

sl_bar1.build()
sl_bar2.build()
sl_bar3.build()
sl_bar4.build()

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
#                   pl_1.diameter/2 *
#                   pl_1.get_asym()[1])

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

qda.get_comp_pos('sublattice_unit_cell', 'sublattice_plunger_1')


# %%

qda.components['sublattice_unit_cell'].cell
qda.components['sublattice_unit_cell'].cell.polygons
qda.components['sublattice_unit_cell'].cell.flatten()
qda.components['sublattice_unit_cell'].cell.polygons


# %%
fog = FanoutGenerator()

(P1, P2, P3, P8, P9, P10) = fog.gates(sl_pl1)
(P4, P5, P6, P7) = fog.gates(sl_pl2)

P1.points
P1.position
P1.layer
P1.add_fanout()

# #%% Generate plungers and barriers

# P1 = sl_pl1.gates[0]
# P2 = sl_pl1.gates[1]
# P3 = sl_pl1.gates[2]
# P4 = sl_pl2.gates[0]
# P5 = sl_pl2.gates[1]
# P6 = sl_pl2.gates[2]
# P7 = sl_pl2.gates[3]
# P8 = sl_pl1.gates[3]
# P9 = sl_pl1.gates[4]
# P10 = sl_pl1.gates[5]


# B1 = sl_bar1[0]
# B2 = sl_bar2[0]
# B3 = sl_bar1[1]
# B4 = sl_bar2[1]
# B5 = sl_bar1[2]
# B6 = sl_bar2[2]
# B7 = sl_bar3[0]
# B8 = sl_bar4[0]
# B9 = sl_bar3[1]
# B10 = sl_bar4[1]
# B11 = sl_bar3[2]
B12 = sl_bar4[2]

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
