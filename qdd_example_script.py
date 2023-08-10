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

bar_45deg = qda_elements.add_barrier('barrier_45deg_rotated')
bar_45deg.width = 30
bar_45deg.length = 70
bar_45deg.layer = 5
bar_45deg.rotate = 1/4*np.pi

bar_135deg = qda_elements.add_copy(bar_45deg, 'barrier_135deg_rotated')
bar_135deg.rotate = 3/4*np.pi

bar_225deg = qda_elements.add_copy(bar_45deg, 'barrier_225deg_rotated')
bar_225deg.rotate = -3/4*np.pi

bar_315deg = qda_elements.add_copy(bar_45deg, 'barrier_315deg_rotated')
bar_315deg.rotate = -1/4*np.pi

bar_45deg.build()
bar_135deg.build()
bar_225deg.build()
bar_315deg.build()

# %% define unit cell

unit_cell.spacing_qd = 200

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
uc_bar_45deg.center = (+qda.spacing_qd_diag/4, qda.spacing_qd_diag/4)
uc_bar_45deg.columns = 2
uc_bar_45deg.spacing = (qda.spacing_qd_diag, qda.spacing_qd_diag)

uc_bar_135deg.component = bar_135deg
uc_bar_135deg.center = (-qda.spacing_qd_diag/4, qda.spacing_qd_diag/4)
uc_bar_135deg.columns = 2
uc_bar_135deg.spacing = (qda.spacing_qd_diag, qda.spacing_qd_diag)

uc_bar_225deg.component = bar_225deg
uc_bar_225deg.center = (-qda.spacing_qd_diag/4, -qda.spacing_qd_diag/4)
uc_bar_225deg.columns = 2
uc_bar_225deg.spacing = (qda.spacing_qd_diag, qda.spacing_qd_diag)

uc_bar_315deg.component = bar_315deg
uc_bar_315deg.center = (+qda.spacing_qd_diag/4, -qda.spacing_qd_diag/4)
uc_bar_315deg.columns = 2
uc_bar_315deg.spacing = (qda.spacing_qd_diag, qda.spacing_qd_diag)

uc_pl_ver.build()
uc_pl_hor.build()

uc_bar_135deg.build()
uc_bar_45deg.build()
uc_bar_315deg.build()
uc_bar_225deg.build()

# %% define sensor

qda_elements.spacing_sep = 60

sensor_top = qda_elements.add_sensor('sensor_top')
sensor_top.source_pos = 'left'
sensor_top.drain_pos = 'right'
sensor_top.sep_pos = 'bottom'
sensor_top.gap_sep = 60
sensor_top.gap_ohmic_pl = 40

sensor_top.plunger.diameter = 160
sensor_top.plunger.layer = 21

sensor_top.barrier_source.width = 30
sensor_top.barrier_source.length = 70
sensor_top.barrier_source.layer = 5

sensor_top.barrier_drain.width = 30
sensor_top.barrier_drain.length = 70
sensor_top.barrier_drain.layer = 5

sensor_top.barrier_sep.width = 50
sensor_top.barrier_sep.length = 60
sensor_top.barrier_sep.layer = 5

sensor_bottom = qda_elements.add_copy(sensor_top, 'sensor_bottom')
sensor_bottom.source_pos = 'right'
sensor_bottom.drain_pos = 'left'
sensor_bottom.sep_pos = 'top'

sensor_right = qda_elements.add_copy(sensor_top, 'sensor_right')
sensor_right.source_pos = 'top'
sensor_right.drain_pos = 'bottom'
sensor_right.sep_pos = 'left'

sensor_left = qda_elements.add_copy(sensor_top, 'sensor_left')
sensor_left.source_pos = 'bottom'
sensor_left.drain_pos = 'top'
sensor_left.sep_pos = 'right'

sensor_top.build_elements()
sensor_top.build()
sensor_bottom.build_elements()
sensor_bottom.build()
sensor_right.build_elements()
sensor_right.build()
sensor_left.build_elements()
sensor_left.build()

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
uc_unitcell.columns = 10
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

# %% Build and save

qda.build()
qda.save_as_gds('qdd_test_design.gds')


# %%
fog = FanoutGenerator()

(P1, P2, P3, P8, P9, P10) = fog.gates(uc_pl_ver)
(P4, P5, P6, P7) = fog.gates(uc_pl_hor)

P1.points
P1.position
P1.layer
P1.add_fanout()

# #%% Generate plungers and barriers

# P1 = uc_pl_ver.gates[0]
# P2 = uc_pl_ver.gates[1]
# P3 = uc_pl_ver.gates[2]
# P4 = uc_pl_hor.gates[0]
# P5 = uc_pl_hor.gates[1]
# P6 = uc_pl_hor.gates[2]
# P7 = uc_pl_hor.gates[3]
# P8 = uc_pl_ver.gates[3]
# P9 = uc_pl_ver.gates[4]
# P10 = uc_pl_ver.gates[5]


# B1 = uc_bar_135deg[0]
# B2 = uc_bar_45deg[0]
# B3 = uc_bar_135deg[1]
# B4 = uc_bar_45deg[1]
# B5 = uc_bar_135deg[2]
# B6 = uc_bar_45deg[2]
# B7 = uc_bar_315deg[0]
# B8 = uc_bar_225deg[0]
# B9 = uc_bar_315deg[1]
# B10 = uc_bar_225deg[1]
# B11 = uc_bar_315deg[2]
B12 = uc_bar_225deg[2]

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
