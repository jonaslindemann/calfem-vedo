# -*- coding: iso-8859-15 -*-

import sys
import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv

import vedo as vd
import vis_vedo as vvd

import beam_solver as bs
import beam_vis as bv

# --- Import model

import dome as cfmodel
#import bridge_beam_bar as cfmodel

if __name__ == "__main__":
    
    solver = bs.BeamBarSolver(cfmodel)
    solver.execute()

    vis = bv.BeamVisualiser(solver)
    vis.draw()
    vis.show_and_wait()






