# -*- coding: iso-8859-15 -*-

import sys
import numpy as np
import calfem.core as cfc

def calc_limits(coords):

    if coords.shape[0]>0:

        max_x = np.max(coords[:,0])
        min_x = np.min(coords[:,0])
        max_y = np.max(coords[:,1])
        min_y = np.min(coords[:,1])
        max_z = np.max(coords[:,2])
        min_z = np.min(coords[:,2])

        return max_x, min_x, max_y, min_y, max_z, min_z
    else:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

def calc_size(coords):

    max_x, min_x, max_y, min_y, max_z, min_z = calc_limits(coords)
    lx = max_x - min_x
    ly = max_y - min_y
    lz = max_z - min_z

    return lx, ly, lz

def calc_beam_displ_limits(a, coords, edof, dofs):

    if edof.shape[0]>0:

        ex, ey, ez = cfc.coord_extract(edof, coords, dofs)
        ed = cfc.extract_eldisp(edof, a)

        coords, topo, node_dofs = convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node=6, ignore_first=False)

        min_displ = 1e300
        max_displ = -1e300

        for el_topo in topo:
            d0 = np.array(a[node_dofs[el_topo[0]][:3]-1]).flatten()
            d1 = np.array(a[node_dofs[el_topo[1]][:3]-1]).flatten()

            l_d0 = np.linalg.norm(d0) 
            l_d1 = np.linalg.norm(d1) 

            if l_d0>max_displ:
                max_displ = l_d0
            if l_d1>max_displ:
                max_displ = l_d1
            if l_d0<min_displ:
                min_displ = l_d0
            if l_d1<min_displ:
                min_displ = l_d1

        return min_displ, max_displ

    else: 

        return 0.0, 0.0
    
def calc_bar_displ_limits(a, coords, edof, dofs):

    if edof.shape[0]>0:

        ex, ey, ez = cfc.coord_extract(edof, coords, dofs)
        ed = cfc.extract_eldisp(edof, a)

        coords, topo, node_dofs = convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node=3, ignore_first=False)

        min_displ = 1e300
        max_displ = -1e300

        for el_topo in topo:
            d0 = np.array(a[node_dofs[el_topo[0]][:3]-1]).flatten()
            d1 = np.array(a[node_dofs[el_topo[1]][:3]-1]).flatten()

            l_d0 = np.linalg.norm(d0)
            l_d1 = np.linalg.norm(d1)

            if l_d0>max_displ:
                max_displ = l_d0
            if l_d1>max_displ:
                max_displ = l_d1
            if l_d0<min_displ:
                min_displ = l_d0
            if l_d1<min_displ:
                min_displ = l_d1

        return min_displ, max_displ

    else:
        
        return 0.0, 0.0
    
def convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node=3, ignore_first=True):
    """
    Routine to convert dof based topology and element coordinates to node based
    topology required for visualisation with VTK and other visualisation frameworks

    :param array edof: element topology [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]
    :param array ex: element x coordinates [nel x n_nodes]
    :param array ey: element y coordinates [nel x n_nodes]
    :param array ez: element z coordinates [nel x n_nodes]
    :param array n_dofs_per_node: number of dofs per node. (default = 3)
    :param boolean ignore_first: ignore first column of edof. (default = True)
    :return array coords: Array of node coordinates. [n_nodes x 3]
    :return array topo: Node topology. [nel x n_nodes]
    :return array node_dofs: Dofs for each node. [n_nodes x n_dofs_per_node]
    """

    node_hash_coords = {}
    node_hash_numbers = {}
    node_hash_dofs = {}
    el_hash_dofs = []

    nel, cols = edof.shape

    if ignore_first:
        tot_dofs = cols-1
    else:
        tot_dofs = cols

    n_nodes = int(tot_dofs / n_dofs_per_node)

    for elx, ely, elz, dofs in zip(ex, ey, ez, edof):

        if ignore_first:
            el_dofs = dofs[1:]
        else:
            el_dofs = dofs

        # 0 1 2  3 4 5  6 7 8  9 12 11 

        el_dof = np.zeros((n_nodes, n_dofs_per_node), dtype=int)
        el_hash_topo = []

        for i in range(n_nodes):
            el_dof[i] = el_dofs[ (i*n_dofs_per_node):((i+1)*n_dofs_per_node) ]
            node_hash_coords[hash(tuple(el_dof[i]))] = [elx[i], ely[i], elz[i]]
            node_hash_numbers[hash(tuple(el_dof[i]))] = -1
            node_hash_dofs[hash(tuple(el_dof[i]))] = el_dof[i]
            el_hash_topo.append(hash(tuple(el_dof[i])))

        el_hash_dofs.append(el_hash_topo)

    coord_count = 0

    coords = []
    node_dofs = []

    for node_hash in node_hash_numbers.keys():
        node_hash_numbers[node_hash] = coord_count
        node_dofs.append(list(node_hash_dofs[node_hash]))
        coord_count +=1

        coords.append(node_hash_coords[node_hash])

    topo = []

    for el_hashes in el_hash_dofs:
        el_topo = []
        for i in range(n_nodes):
            el_topo.append(node_hash_numbers[el_hashes[i]])

        topo.append(el_topo)

    return np.array(coords), np.array(topo), np.array(node_dofs)
