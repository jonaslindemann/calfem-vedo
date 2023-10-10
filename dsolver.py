# -*- coding: iso-8859-15 -*-
'''
ObjectiveFrame generated CALFEM Python code.
'''

import numpy as np
import calfem.core as cfc

import vedo as vd
import vis_vedo as vvd

# --- Import model

from bridge_beam_bar import *

# --- Extract element coordinates

if edof_beams.shape[0]>0:
    ex_bm, ey_bm, ez_bm = cfc.coord_extract(edof_beams, coords_beams, dofs_beams)

if edof_bars.shape[0]>0:
    ex_br, ey_br, ez_br = cfc.coord_extract(edof_bars, coords_bars, dofs_bars)

# --- Setup global matrices

max_dof_beams = -1
max_dof_bars = -1

if edof_beams.shape[0]>0:
    max_dof_beams = np.max(edof_beams)

if edof_bars.shape[0]>0:
    max_dof_bars = np.max(edof_bars)

n_dofs = max(max_dof_bars, max_dof_beams)

n_beams = edof_beams.shape[0]
n_bars = edof_bars.shape[0]

K = np.zeros((n_dofs, n_dofs))
f = np.zeros((n_dofs, 1))

# --- Assemble beams

if n_beams>0:
    for elx, ely, elz, dofs, epi, eo in zip(ex_bm, ey_bm, ez_bm, edof_beams, beam_mat_idx, beam_orientation):
        Ke = cfc.beam3e(elx, ely, elz, eo, ep[epi])
        K = cfc.assem(dofs, K, Ke)

# --- Assemble bars

if n_bars>0:
    for elx, ely, elz, dofs, epi in zip(ex_br, ey_br, ez_br, edof_bars, bar_mat_idx):
        E = ep[epi][0]
        A = ep[epi][2]
        Ke = cfc.bar3e(elx, ely, elz, [E, A])
        K = cfc.assem(dofs, K, Ke)

# --- Setting up loads and boundary conds

if False:
    bc = []
    bc_vals = []

    for ni, bi in node_bc_idx:
        bc_prescr_dofs = node_bc[bi]
        bc_values = node_bc_val[bi]
        bc_dofs = dofs_beams[ni-1]
        for i, p_dof in enumerate(bc_prescr_dofs):
            if p_dof == 1:
                bc.append(bc_dofs[i])
                bc_vals.append(bc_values[i])

    for ni, li in node_load_idx:
        load_values = node_load[li]
        load_dofs = dofs_beams[ni-1][0:3]
        f[load_dofs-1, 0] += load_values

    print(f.sum())

    bc = np.array(bc)
    bc_vals = np.array(bc_vals)

    # --- Solve equation system

    a, r = cfc.solveq(K, f, bc, bc_vals)

    # --- Extract element forces

    ed_beams = cfc.extract_ed(edof_beams, a)

    beam_forces = []

    for elx, ely, elz, ed, epi, eo in zip(ex_bm, ey_bm, ez_bm, ed_beams, beam_mat_idx, beam_orientation):
        es, edi, eci = cfc.beam3s(elx, ely, elz, eo, ep[epi], ed, nep=10)
        beam_forces.append([es, edi, eci])

    print(r.sum())

if n_bars>0:
    vvd.draw_bars(edof_bars, ex_br, ey_br, ez_br)

if n_beams>0:
    vvd.draw_beams(edof_beams, ex_bm, ey_bm, ez_bm)



vvd.show_and_wait()


