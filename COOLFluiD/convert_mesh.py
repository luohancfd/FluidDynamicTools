# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 10:51:16 2018
A script to convert Gmesh to CFmesh

@author: Han
"""
import os
import meshio
import numpy as np
from collections import OrderedDict


num_nodes_per_cell = {
    "vertex": 1,
    "line": 2,
    "triangle": 3,
    "quad": 4,
    "quad8": 8,
    "tetra": 4,
    "hexahedron": 8,
    "hexahedron20": 20,
    "wedge": 6,
    "pyramid": 5,
    #
    "line3": 3,
    "triangle6": 6,
    "quad9": 9,
    "tetra10": 10,
    "hexahedron27": 27,
    "wedge18": 18,
    "pyramid14": 14,
    #
    "line4": 4,
    "triangle10": 10,
    "quad16": 16,
    "tetra20": 20,
    "wedge40": 40,
    "hexahedron64": 64,
    #
    "line5": 5,
    "triangle15": 15,
    "quad25": 25,
    "tetra35": 35,
    "wedge75": 75,
    "hexahedron125": 125,
    #
    "line6": 6,
    "triangle21": 21,
    "quad36": 36,
    "tetra56": 56,
    "wedge126": 126,
    "hexahedron216": 216,
    #
    "line7": 7,
    "triangle28": 28,
    "quad49": 49,
    "tetra84": 84,
    "wedge196": 196,
    "hexahedron343": 343,
    #
    "line8": 8,
    "triangle36": 36,
    "quad64": 64,
    "tetra120": 120,
    "wedge288": 288,
    "hexahedron512": 512,
    #
    "line9": 9,
    "triangle45": 45,
    "quad81": 81,
    "tetra165": 165,
    "wedge405": 405,
    "hexahedron729": 729,
    #
    "line10": 10,
    "triangle55": 55,
    "quad100": 100,
    "tetra220": 220,
    "wedge550": 550,
    "hexahedron1000": 1000,
    #
    "line11": 11,
    "triangle66": 66,
    "quad121": 121,
    "tetra286": 286,
}


class Node:
    def __init__(self, id):
        self.id = id
        self.cellID = set([])
        self.bcID = set([])

    def setCell(self, cellID):
        self.cellID.add(cellID)

    def setBc(self, bcID):
        self.bcID.add(bcID)

    def __repr__(self):
        return 'ID: %d\nBC_ID: %s\nCELL_ID: %s' % (self.id, str(self.bcID), str(self.cellID))

    def copy(self, b):
        self.cellID = b.cellID.copy()
        self.bcID = b.bcID.copy()


class Nodes:
    '''
    A collection of nodes
    '''

    def __init__(self, Nnode):
        self.Nodes = [Node(i) for i in range(Nnode)]

    def __getitem__(self, a):
        return self.Nodes[a]

    def __setitem__(self, a, b):
        self.Nodes[a].copy(b)

    def __repr__(self):
        return 'Number of node: %d' % (len(self.Nodes))

    def __len__(self):
        return len(Nodes)

    def setInfo(self, data, settype, ID=None):
        '''
        set the nodes cellID/bcID

        data: a numpy matrix, each row is the id of nodes
        settype: 'cell' or 'bc'
        ID(opt): a list, id of the cell/bc
        '''
        if not ID:
            if type(data) == np.ndarray:
                ID = range(data.shape[0])
            else:
                ID = range(len(data))

        if settype == 'cell':
            for i, j in zip(data, ID):
                for k in i:
                    self.Nodes[k].setCell(j)
        elif settype == 'bc':
            for i, j in zip(data, ID):
                for k in i:
                    self.Nodes[k].setBc(j)

    def exportInfo(self, settype, index=None):
        result = []
        if not index:
            index = range(len(self.Nodes))

        if settype == 'cell':
            for i in index:
                result.append(sorted(list(self.Nodes[i].cellID)))
        elif settype == 'bc':
             for i in index:
                result.append(sorted(list(self.Nodes[i].bcID)))
        return result


def nparray2string(val):
    '''
    convert a numpy array to string
    '''
    return ' ' + np.array2string(val, threshold=np.nan,
                                 max_line_width=np.inf,
                                 formatter={'int': lambda x: '%d' % (x)}).replace('[', '').replace(']', '') + '\n'


class Cells:
    '''
    storage of cell nodes
    '''

    def __init__(self, data, dim=2, mesh_type='cc'):
        __cell_type__ = [
            {1: 'point'},  # 0d
            {2: 'line'},   # 1d
            {3: 'triangle', 4: 'quad',    5: 'pyramid'},  # 2d
            {4: 'tetra',    5: 'pyramid', 6: 'wedge'},  # 3d
        ]
        self.dim = dim   # dimension of the grid
        self.mesh_type = mesh_type
        if type(data) == np.ndarray:
            self.data = [data]
            self.cell_type = [__cell_type__[dim - 1][len(data[0])]]
        elif type(data) == dict:
            self.data = list(data.values())
            self.cell_type = list(data.keys())
        self.num_cell_type = [i.shape[0] for i in self.data]  # number of cells for each type
        self.num_cell = sum(self.num_cell_type)  # total number of cells
        self.cell_dim = list(range(len(self.num_cell_type)))
        for i, k in enumerate(self.cell_type):
            for l, j in enumerate(__cell_type__):
                if k in list(j.values()):
                    self.cell_dim[i] = l
        # put elements in front of bc
        sort_index = [i[0] for i in sorted(
            enumerate(self.cell_dim), key=lambda x: x[1], reverse=True)]
        self.data = [self.data[i] for i in sort_index]
        self.num_cell_type = [self.num_cell_type[i] for i in sort_index]
        self.cell_type = [self.cell_type[i] for i in sort_index]
        self.cell_dim = [self.cell_dim[i] for i in sort_index]

        # assign id
        self.cellID = []
        self.stride = [0] + [len(j) for i, j in enumerate(self.data[:-1])]
        self.end = []
        for j in range(len(self.stride)):
            if j > 0:
                self.stride[j] += self.stride[j - 1]
        self.end = [i for i in self.stride[1:]]
        self.end = self.end + [self.num_cell]

        # create nodes, first calculate number of nodes
        self.num_node = 0
        for i in self.data:
            j = np.max(i)
            if j > self.num_node:
                self.num_node = j
        self.num_node += 1
        self.nodes = Nodes(self.num_node)

        # set cell info and bc info for the nodes
        self.num_bc = 0
        for i, (nodes, cell_type, cell_dim) in enumerate(zip(self.data, self.cell_type, self.cell_dim)):
            if cell_dim > self.dim - 1:
                self.nodes.setInfo(nodes, 'cell', ID=list(
                    range(self.stride[i], self.end[i])))
            else:
                self.nodes.setInfo(nodes, 'bc', ID=range(
                    self.stride[i], self.end[i]))

        # set stateID
        self.stateID = []
        if mesh_type == 'cc':
            for imesh, (s, e, d) in enumerate(zip(self.stride, self.end, self.cell_dim)):
                stateID = np.arange(s, e)  # surface cell
                if d < self.dim:
                    # bc cell
                    for icell, nodes in enumerate(self.data[imesh]):
                        w = set.intersection(
                            *[self.nodes[inode].cellID for inode in nodes])
                        if len(w) == 1:
                            stateID[icell] = w.pop()
                        else:
                            raise ValueError(
                                'More than one cell have the same nodes')
                self.stateID.append(stateID)
        elif mesh_type == 'vc':
            for i in self.data:
                self.stateID.append(i)

        # a list containing index of bc, sorted to ensure correct rank
        # use readBC(self, cell_data, field_data) to set the data
        self.cell_phys_type = [None for i in self.data]

    def __getnum__(self, id):
        # get two index to search node from data
        l = len(self.stride) - 1
        for i, j in enumerate(self.stride):
            if id < j:
                l = i - 1
                break
        l2 = id - self.stride[l]
        return l, l2

    def __getitem__(self, id):
        l, l2 = self.__getnum__(id)
        return self.data[l][l2]

    def __setitem__(self, id, item):
        l, l2 = self.__getnum__(id)
        assert len(self.data[l][l2]) != len(item)
        self.data[l][l2] = item

    def exportCFmeshCell(self, index=None):
        result = []
        if not index:
            index = list(range(len(self.data)))

        for i in index:
            if self.cell_dim[i] >= self.dim:  # element cell
                if self.stateID[i].ndim == 1:
                    w = np.vstack((self.data[i].T, self.stateID[i])).T
                else:
                    w = np.vstack((self.data[i].T, self.stateID[i].T)).T
                result.append(nparray2string(w))
            # not test yet
            else:
                col1 = np.ones(
                    (self.num_cell_type[i]), dtype=np.int)*self.data[i].shape[1]
                if self.stateID[i].ndim == 1:
                    col2 = np.ones((self.num_cell_type[i]), dtype=np.int)
                else:
                    col2 = np.ones(
                        (self.num_cell_type[i]), dtype=np.int)*self.stateID[i].shape[1]
                w = np.vstack(
                    (col1, col2, self.data[i].T, self.stateID[i].T)).T
                result.append(nparray2string(w))
        return result

    def exportCFmeshBC(self, index=None):
        '''
        Use to export Line BC, for surface BC, use exportCFmeshCell

        index: a list containing the index of BCs to export
        '''

        assert self.bc is not None
        if not index:
            index = range(len(self.bc))

        result = ['' for i in index]
        for ibc in index:
            bc = self.bc[ibc]
            if not bc:
                raise ValueError('No bc info is found for %d' % (ibc,))

            nTRs = len(bc)
            header = '!TRS_NAME %s\n!NB_TRs %d\n!NB_GEOM_ENTS '%(self.bc_type_name[self.bc_type[ibc]], nTRs)
            w2 = ['' for i in range(nTRs)]
            for ibc2, bc2 in enumerate(bc):
                if ibc2 < len(bc)-1:
                    header += '%d ' % (len(bc2))
                else:
                    header += '%d\n!GEOM_TYPE Face\n!LIST_GEOM_ENT\n' % (len(bc2))
                # bc2 store the id of lines
                if len(bc2) < 1:
                    raise ValueError('No enough points in %d->%d'(ibc, ibc2))
                l, l2 = self.__getnum__(bc2[0])

                if self.stateID[l].ndim == 1:
                    col12 = np.array([2, 1])
                else:
                    col12 = np.array([2, 2])

                w3 = np.hstack((col12, self.data[l][l2].flatten(), self.stateID[l][l2].flatten()))

                for iline, line_id in enumerate(bc2[1:]):
                    l, l2 = self.__getnum__(line_id)
                    w4 = np.hstack((col12, self.data[l][l2].flatten(), self.stateID[l][l2].flatten()))
                    w3 = np.vstack((w3, w4))

                # convert the matrix to string
                w2[ibc2] = nparray2string(w3)

            result[ibc] = ''.join([header] + w2)
        return result

    def __repr__(self):
        return 'Cell type: ' + str(self.cell_type)+'\nNumber of cells: ' + str([i.shape[0] for i in self.data]) + '\nNumber of nodes: %d' % (self.num_node)

    def printCellInfo(self, inode):
        l, l2 = self.__getnum__(inode)
        nodes = 'Nodes: ' + str(self.data[l][l2])
        cellID = 'CellID: %d' % (inode,)
        stateID = 'StateID: ' + str(self.stateID[l][l2])
        print(cellID + '\n' + stateID + '\n' + nodes + '\n')

    def readBC(self, cell_data, field_data):
        '''
        Load bc/cell type number

        cell_data: mesh.cell_data
        '''
        for i, dat in enumerate(cell_data.values()):
            for key, value in dat.items():
                if 'physical' in key:
                    self.cell_phys_type[i] = value
                    break
        self.bc_type = set([])
        for cell_phys_type, cell_dim in zip(self.cell_phys_type, self.cell_dim):
            if cell_dim < self.dim:
                self.bc_type = self.bc_type.union(set(cell_phys_type))
        self.nbc_type = len(self.bc_type)
        self.bc_type = list(self.bc_type)
        self.bc_type_name = {i[0]: j for j, i in field_data.items()}

        self.bc = [[] for i in self.bc_type]
        for imesh, nodes in enumerate(self.data):
            if self.cell_dim[imesh] == 1:  # only sort line
                for ibc, bc_type in enumerate(self.bc_type):
                    bc_id = [self.stride[imesh]+i for i,
                        j in enumerate(self.cell_phys_type[imesh]) if j == bc_type]

                    line_list = []
                    for i, iline in enumerate(bc_id):
                        l, l2 = self.__getnum__(iline)
                        line_list.append(
                            Line(iline, [self.nodes[j] for j in self.data[l][l2]]))

                    possible_start = []
                    for i, line in enumerate(line_list):
                        if line.p1 in bc_id:
                            line.p1 = line_list[bc_id.index(line.p1)]
                        else:
                            line.p1 = None
                            possible_start.append(line.id)

                        if line.p2 in bc_id:
                            line.p2 = line_list[bc_id.index(line.p2)]
                        else:
                            line.p2 = None
                            possible_start.append(line.id)

                    # ok, now let's find the one with None, which should be the start or end of line
                    for istart in possible_start:
                        icon = False
                        if self.bc[ibc]:
                            for iibc, w in enumerate(self.bc[ibc]):
                                if istart in w:
                                    icon = True
                                    break
                        if icon:
                            continue
                        else:
                            bc = []
                            i = line_list[bc_id.index(istart)]
                            while i.p1 or i.p2:
                                bc.append(i.id)
                                inext = None
                                if i.p1:
                                    if i.p1.id not in bc:
                                        inext = i.p1
                                if i.p2:
                                    if i.p2.id not in bc:
                                        inext = i.p2
                                if inext:
                                    i = inext
                                else:
                                    break
                            self.bc[ibc].append(bc)


class Line:
    def __init__(self, lineid, nodes):
        self.id = lineid  # id of this line
        # put nodes in desired order
        node_sort = sorted(nodes,  key=lambda x: -x.id)
        ids = [i.id for i in node_sort]

        self.nodes = node_sort
        self.p1 = self.nodes[0].bcID.copy()
        self.p1.remove(lineid)
        self.p2 = self.nodes[1].bcID.copy()
        self.p2.remove(lineid)

        if self.p1:
            self.p1 = list(self.p1)[0]
        else:
            self.p1 = None
        if self.p2:
            self.p2 = list(self.p2)[0]
        else:
            self.p2 = None

    def __repr__(self):
        return str(self.connectInfo())

    def connectInfo(self):
        if not self.p1:
            p1 = None
        elif type(self.p1) == int:
            p1 = self.p1
        else:
            p1 = self.p1.id
        if not self.p2:
            p2 = None
        elif type(self.p2) == int:
            p2 = self.p2
        else:
            p2 = self.p2.id
        return [p1, p2]


def CFmesh_write(mesh, cfmesh_name, nb_dim=2, nb_eq=4, cell_center=True):
    '''
    export CFmesh file

    nb_dim: number of spatial dimension
    nb_eq:  number of equations
    '''
    # write header info
    header = OrderedDict()
    header['COOLFLUID_VERSION'] = '2013.9'
   # header['COOLFLUID_SVNVERSION'] = '1387'
    header['CFMESH_FORMAT_VERSION'] = '1.3'
    header['nb_dim'] = nb_dim
    header['nb_eq'] = nb_eq
    nb_nodes = mesh.points.shape[0]
    nb_elem = 0
    nb_elem_types = 0
    if nb_dim == 2:
        for elem, nodes in mesh.cells.items():
            if 'line' not in elem.lower() and num_nodes_per_cell[elem.lower()] > 2:
                nb_elem += len(nodes)
                nb_elem_types += 1
    else:
        raise ValueError('Does not support nb_dim != 2 yet')

    header['nb_nodes'] = '%d 0' % (nb_nodes,)
    if not cell_center:
        header['nb_states'] = '%d 0' % (nb_nodes,)
    else:
        header['nb_states'] = '%d 0' % (nb_elem,)
    header['nb_elem'] = nb_elem
    header['nb_elem_types'] = nb_elem_types
    header['geom_polyorder'] = 1
    if cell_center:
        header['sol_polyorder'] = 0
    else:
        header['sol_polyorder'] = 1    
    if nb_dim == 2:
        header['elem_types'] = [i.capitalize()
                                             for i in mesh.cells.keys() if 'line' not in i.lower()]
    else:
        raise ValueError('Does not support nb_dim != 2 yet')
    header['nb_elem_per_type'] = [mesh.cells[i.lower()].shape[0]
                                                     for i in header['elem_types']]
    header['nb_nodes_per_type'] = [num_nodes_per_cell[i.lower()]
                                                              for i in header['elem_types']]
    if cell_center:
        header['nb_states_per_type'] = [
            1 for i in range(len(header['elem_types']))]
    else:
        header['nb_states_per_type'] = [num_nodes_per_cell[i.lower()]
                                                                   for i in header['elem_types']]

    if cell_center:
        cell_info = Cells(mesh.cells, dim=2, mesh_type='cc')
    else:
        cell_info = Cells(mesh.cells, dim=2, mesh_type='vc')
    cell_info.readBC(mesh.cell_data, mesh.field_data)

    with open(cfmesh_name, 'w', encoding='utf-8') as f:
        for key, value in header.items():
            if type(value) == str:
                f.write('!%s %s\n' % (key.upper(), value))
            elif type(value) == int:
                f.write('!%s %d\n' % (key.upper(), value))
            elif type(value) == list:
                f.write('!%s %s\n' %
                        (key.upper(), ','.join([str(i) for i in value])))

        # ======== write cell info ====================
        f.write('!LIST_ELEM\n')
        for i, icell_dim in enumerate(cell_info.cell_dim):
            if icell_dim >= cell_info.dim:
                f.write(cell_info.exportCFmeshCell([i])[0])

        # ======== write boundary conditions ==========
        f.write('!NB_TRSs %d\n' % (len(cell_info.bc_type),))
        w = cell_info.exportCFmeshBC()
        for i in w:
            f.write(i)

        # write coordinates of node
        f.write('!LIST_NODE\n')
        f.write(' '+np.array2string(mesh.points[:,:nb_dim],threshold=np.nan,max_line_width=np.inf).
                        replace('[','').replace(']','')+'\n')

        # write state
        f.write('!LIST_STATE 0\n')
        f.write('!END')

def StripGmshComment(gmesh_file):
    gmesh_file2 = 'tmp.msh'
    with open(gmesh_file2, 'w', encoding='utf-8') as f1:
        with open(gmesh_file, 'r', encoding='utf-8') as f2:
            iwrite = True
            for line in f2:
                if '$Comments' in line:
                    iwrite = False
                elif '$EndComments' in line:
                    iwrite = True
                elif iwrite:
                    f1.write(line)
    os.remove(gmesh_file)
    os.rename(gmesh_file2, gmesh_file)

def CRLF2LF(filename):
    with open(filename, "rb") as f:
        data = f.read()
    if b'\0' in data:
        raise ValueError(filename + "is Binary!")
        return
    newdata = data.replace(b"\r\n", b"\n")
    if newdata != data:
        with open(filename, "wb") as f:
            f.write(newdata)    
        
def StripLeadingSpace(cfmesh_name):
    with open(cfmesh_name, 'r') as f:
        content = f.readlines()
    with open(cfmesh_name, 'w') as f:
        for line in content:
            f.write(line.lstrip())

if (__name__ == '__main__'):
    gmesh_file = 'ramp.msh'
    StripGmshComment(gmesh_file)
    mesh = meshio.read(gmesh_file)
    cfmesh_name='ramp.CFmesh'
    CFmesh_write(mesh, cfmesh_name, nb_dim=2, nb_eq=4, cell_center=True)
    StripLeadingSpace(cfmesh_name)


