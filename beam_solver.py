# -*- coding: iso-8859-15 -*-

import sys
import numpy as np
import calfem.core as cfc
import beam_utils as cfu

class BeamBarSolver:
    def __init__(self, cfmodel):
        self.cfmodel = cfmodel

        self.edof_beams = cfmodel.edof_beams
        self.edof_bars = cfmodel.edof_bars
        self.coords_beams = cfmodel.coords_beams
        self.coords_bars = cfmodel.coords_bars
        self.dofs_bars = cfmodel.dofs_bars
        self.dofs_beams = cfmodel.dofs_beams
        self.beam_orientation = cfmodel.beam_orientation
        self.beam_mat_idx = cfmodel.beam_mat_idx
        self.bar_mat_idx = cfmodel.bar_mat_idx
        self.ep = cfmodel.ep
        self.node_bc = cfmodel.node_bc
        self.node_bc_val = cfmodel.node_bc_val
        self.node_bc_idx = cfmodel.node_bc_idx
        self.node_load = cfmodel.node_load
        self.node_load_idx = cfmodel.node_load_idx

        self.magnfac_percent = 0.1

    def setup(self):

        # --- Extract element coordinates

        if self.edof_beams.shape[0]>0:
            self.ex_bm, self.ey_bm, self.ez_bm = cfc.coord_extract(self.edof_beams, self.coords_beams, self.dofs_beams)

        if self.edof_bars.shape[0]>0:
            self.ex_br, self.ey_br, self.ez_br = cfc.coord_extract(self.edof_bars, self.coords_bars, self.dofs_bars)

        # --- Setup global matrices

        max_dof_beams = -1
        max_dof_bars = -1

        if self.edof_beams.shape[0]>0:
            max_dof_beams = np.max(self.edof_beams)

        if self.edof_bars.shape[0]>0:
            max_dof_bars = np.max(self.edof_bars)

        self.n_dofs = max(max_dof_bars, max_dof_beams)

        self.n_beams = self.edof_beams.shape[0]
        self.n_bars = self.edof_bars.shape[0]

        self.K = np.zeros((self.n_dofs, self.n_dofs))
        self.f = np.zeros((self.n_dofs, 1))

        # --- Assemble beams

        if self.n_beams>0:
            for elx, ely, elz, dofs, epi, eo in zip(self.ex_bm, self.ey_bm, self.ez_bm, self.edof_beams, self.beam_mat_idx, self.beam_orientation):
                Ke = cfc.beam3e(elx, ely, elz, eo, self.ep[epi])
                self.K = cfc.assem(dofs, self.K, Ke)

        # --- Assemble bars

        if self.n_bars>0:
            for elx, ely, elz, dofs, epi in zip(self.ex_br, self.ey_br, self.ez_br, self.edof_bars, self.bar_mat_idx):
                E = self.ep[epi][0]
                A = self.ep[epi][2]
                Ke = cfc.bar3e(elx, ely, elz, [E, A])
                self.K = cfc.assem(dofs, self.K, Ke)

        # --- Setting up loads and boundary conds

        bc = []
        bc_vals = []

        for ni, bi in self.node_bc_idx:
            bc_prescr_dofs = self.node_bc[bi]
            bc_values = self.node_bc_val[bi]
            bc_dofs = self.dofs_beams[ni-1]
            for i, p_dof in enumerate(bc_prescr_dofs):
                if p_dof == 1:
                    bc.append(bc_dofs[i])
                    bc_vals.append(bc_values[i])

        for ni, li in self.node_load_idx:
            load_values = self.node_load[li]
            load_dofs = self.dofs_beams[ni-1][0:3]
            self.f[load_dofs-1, 0] += load_values

        self.bc = np.array(bc)
        self.bc_vals = np.array(bc_vals)

    def solve(self):                     

        # --- Solve equation system

        self.a, self.r = cfc.solveq(self.K, self.f, self.bc, self.bc_vals)

    def extract_forces(self):

        # --- Extract element forces

        if self.n_beams>0:
            self.ed_beams = cfc.extract_ed(self.edof_beams, self.a)

        if self.n_bars>0:
            self.ed_bars = cfc.extract_ed(self.edof_bars, self.a)

        self.beam_forces = []

        if self.n_beams>0:
            for elx, ely, elz, ed, epi, eo in zip(self.ex_bm, self.ey_bm, self.ez_bm, self.ed_beams, self.beam_mat_idx, self.beam_orientation):
                es, edi, eci = cfc.beam3s(elx, ely, elz, eo, self.ep[epi], ed, nep=10)
                self.beam_forces.append([es, edi, eci])

        self.bar_forces = []

        if self.n_bars>0:
            for elx, ely, elz, ed, epi in zip(self.ex_br, self.ey_br, self.ez_br, self.ed_bars, self.bar_mat_idx):
                es, edi, eci = cfc.bar3s(elx, ely, elz, self.ep[epi][:2], ed, nep=10)
                self.bar_forces.append([es, edi, eci])


        min_displ_bm, max_displ_bm = cfu.calc_beam_displ_limits(self.a, self.coords_beams, self.edof_beams, self.dofs_beams)
        min_displ_br, max_displ_br = cfu.calc_bar_displ_limits(self.a, self.coords_bars, self.edof_bars, self.dofs_bars)

        self.max_displ = max(max_displ_bm, max_displ_br)
        self.min_displ = min(min_displ_bm, min_displ_br)

        lx_bm, ly_bm, lz_bm = cfu.calc_size(self.coords_beams)
        lx_br, ly_br, lz_br = cfu.calc_size(self.coords_bars)

        self.lx = max(lx_bm, lx_br)
        self.ly = max(ly_bm, ly_br)
        self.lz = max(lz_bm, lz_br)

        print("max/min displ ", self.max_displ, self.min_displ)
        print("model size ", self.lx, self.ly, self.lz)

        self.magnfac = max(self.lx, self.ly, self.lz)*self.magnfac_percent/self.max_displ

    def execute(self):
        self.setup()
        self.solve()
        self.extract_forces()
