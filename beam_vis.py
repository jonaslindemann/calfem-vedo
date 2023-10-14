# -*- coding: iso-8859-15 -*-

import sys
import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv

import vedo as vd
import vis_vedo as vvd

import beam_solver as bs

class BeamVisualiser:
    def __init__(self, solver):
        self.solver = solver

    def draw_geometry(self):

        if self.solver.n_bars>0:
            vvd.set_color_and_alpha("yellow", 1.0)
            vvd.draw_bars(self.solver.edof_bars, self.solver.ex_br, self.solver.ey_br, self.solver.ez_br)

        if self.solver.n_beams>0:
            vvd.set_color_and_alpha("yellow", 1.0)
            vvd.draw_beams(self.solver.edof_beams, self.solver.ex_bm, self.solver.ey_bm, self.solver.ez_bm)

    def draw_displacements(self):

        if self.solver.n_bars>0:
            vvd.set_color_and_alpha("gray", 0.7)
            vvd.draw_bar_displacements(self.solver.a, self.solver.coords_bars, self.solver.edof_bars, self.solver.dofs_bars, self.solver.magnfac)

        if self.solver.n_beams>0:
            vvd.set_color_and_alpha("gray", 0.7)
            vvd.draw_beam_displacements(self.solver.a, self.solver.coords_beams, self.solver.edof_beams, self.solver.dofs_beams, self.solver.magnfac)

    def draw(self):
        self.draw_geometry()
        self.draw_displacements()

    def show_and_wait(self):
        vvd.show_and_wait()

