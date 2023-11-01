# -*- coding: iso-8859-15 -*-

import beam_solver as bs
import beam_vis as bv

# --- Import model

import dome as model
#import bridge_beam as model
#import bridge_bar as model
#import bridge_beam_bar as model

if __name__ == "__main__":
    
    solver = bs.BeamBarSolver(model)
    solver.execute()

    vis = bv.BeamVisualiser(solver)
    vis.draw()
    vis.show_and_wait()






