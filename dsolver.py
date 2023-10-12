# -*- coding: iso-8859-15 -*-

import sys
import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv

import vedo as vd
import vis_vedo as vvd

import beam_solver as bs

# --- Import model

#import dome as cfmodel
import bridge_beam_bar as cfmodel

if __name__ == "__main__":
    
    solver = bs.BeamBarSolver(cfmodel)
    solver.execute()


    # if n_bars>0:
    #     vvd.set_color_and_alpha("yellow", 1.0)
    #     vvd.draw_bars(edof_bars, ex_br, ey_br, ez_br)

    #     vvd.set_color_and_alpha("gray", 0.7)
    #     vvd.draw_bar_displacements(a, coords_bars, edof_bars, dofs_bars, magnfac=magnfac)

    # if n_beams>0:
    #     vvd.set_color_and_alpha("yellow", 1.0)
    #     vvd.draw_beams(edof_beams, ex_bm, ey_bm, ez_bm)

    #     vvd.set_color_and_alpha("gray", 0.7)
    #     vvd.draw_beam_displacements(a, coords_beams, edof_beams, dofs_beams, magnfac=magnfac)

    # vvd.show_and_wait()


