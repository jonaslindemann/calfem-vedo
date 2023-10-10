# -*- coding: iso-8859-15 -*-

import numpy as np
import calfem.core as cfc
import vedo as vd

class VedoCanvas:
    """
    Class for handling the VedoMainWindow class and making sure there is only one instance of it.
    """
    __instance = None
    @staticmethod
    def instance():
        """ Static access method. """
        if VedoCanvas.__instance == None: VedoCanvas()
        
        return VedoCanvas.__instance

    def __init__(self):
        """ Virtually private constructor. """
        if VedoCanvas.__instance is None:
            VedoCanvas.__instance = self
            self.__entities = []

    def add(self, vedo_object):
        self.__entities.append(vedo_object)

    @property
    def entities(self):
        return self.__entities
    
    def save(self, filename):
        print(self.__entities)
        vd.write(self.__entities, filename)

def add_canvas(vedo_object):
    VedoCanvas.instance().add(vedo_object)

def save(filename):
    VedoCanvas.instance().save(filename)


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

    print("cols    =", tot_dofs)
    print("nel     =", nel)
    print("n_nodes =", n_nodes)

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

    return coords, topo, node_dofs


def draw_beams(edof, ex, ey, ez):
    coords, topo, node_dofs = convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node=6, ignore_first=False)

    lines = []

    for el_topo in topo:
        l = vd.Line(coords[el_topo[0]], coords[el_topo[1]])
        l.linewidth(2)
        t = vd.Tube([coords[el_topo[0]], coords[el_topo[1]]], r=0.025, c="yellow")
        add_canvas(t)

def draw_bars(edof, ex, ey, ez):
    coords, topo, node_dofs = convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node=3, ignore_first=False)

    lines = []

    for el_topo in topo:
        l = vd.Line(coords[el_topo[0]], coords[el_topo[1]])
        l.linewidth(2)
        t = vd.Tube([coords[el_topo[0]], coords[el_topo[1]]], r=0.025, c="orange")
        add_canvas(t)

def show_and_wait():
    plt = vd.Plotter(axes=1, bg2='lightblue')
    plt.look_at(plane="xy")
    plt.show(VedoCanvas.instance().entities, __doc__, viewup='z')
    plt.close()

