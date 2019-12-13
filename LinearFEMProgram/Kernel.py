"""
Created on Fri Nov 01 2019
Main Programme for linear elastic analysis under concentrated loads;
reads mesh, bcds and loads from the part definition of 
Abaqus input file (only a single part is supported);
if the input file is created by Abaqus, then both nset and elset are supported;
if the input file is created by Salome, then only nset is supported (elset 
requires renumbering of elements due to unused boundary elements generated 
by Salome)
This module is called by the job file (e.g., Job-1.py)
@author: Boyang CHEN TU Delft
"""
import numpy as np
from . import Parameters as param
from .Preprocessing import read_abaqus_parts as read_abaqus
from .Elements import frame2d2elem, tri2d3elem
from .Solvers import assembler, direct_solver_edu as solver

# =================================================
# Define classes for objects in kernel program
# =================================================
class dimension_data:
    def __init__(self, NDIM, NST, NDOF_NODE, ELEM_TYPES):
        self.NDIM = NDIM
        self.NST  = NST
        self.NDOF_NODE = NDOF_NODE
        self.ELEM_TYPES = ELEM_TYPES

class nset_data:
    def __init__(self, settype, nodalDofs, dofValues):
        self.settype = settype # either 'bcd' or 'cload'
        self.nodalDofs = nodalDofs # indices of affected nodal DoFs
        self.dofValues = dofValues # values of imposed disp or cload
        
class dload_function:
    def __init__(self, expression, coord_system='global', order=0):
        self.expression = expression
        self.coord_system = coord_system
        self.order = order

class dload_data:
    def __init__(self, elset, function):
        self.elset = elset
        self.function = function

class elem_list:
    def __init__(self, eltype, elems=None):
        self.eltype = eltype
        self.elems  = elems




def kernel_program(inputfile, dimData, Materials, dict_elset_matID, \
                   dict_nset_data, dload_functions=[], dict_elset_dloadID={}):
    ###############################################################################
    # Preprocessing 
    ###############################################################################
    # Read data from Abaqus input file and form abaqus parts
    parts = read_abaqus.read_parts_from_inputfile(inputfile)
    
    # check if there is only one part
    # in the future, consider making a loop over all parts
    if(not len(parts)==1):
        raise ValueError('Only a single part is supported!')
    
    # verification of dimensional parameters before proceeding
    verify_dimensional_parameters(parts[0], dimData)
        
    # form lists of nodes and elems
    nodes = form_nodes(parts[0])
    elem_lists = form_elem_lists(parts[0], dimData.NDOF_NODE, dimData.ELEM_TYPES, dict_elset_matID)
    
    # form lists of bcds and cloads
    [bcd_dofs, bcd_values, cload_dofs, cload_values] = \
    form_bcds_cloads(parts[0], dict_nset_data, dimData.NDOF_NODE)
    
    # form lists of elset for distributed loads
    list_dload_data = form_list_dload_data(parts[0], dict_elset_dloadID, dload_functions)

    
    ###############################################################################
    # Assembler 
    # obtain the full stiffness matrix K and external distributed force vector f
    ###############################################################################
    [K, f] = assembler(nodes, elem_lists, dimData.NDOF_NODE, Materials, list_dload_data)
    
    
    ###############################################################################
    # Solver
    # modify the stiffness matrix and force vector based on applied bcds and cloads
    # obtain dof vector a and reaction force vector RF, both size ndof by 1
    ###############################################################################    
    [a, RF] = solver(K, f, bcd_dofs, bcd_values, cload_dofs, cload_values)
    
    
    return [parts, nodes, elem_lists, f, a, RF]





def verify_dimensional_parameters(part, dimData):
    
    for elem_group in part.elem_groups:
        # check if ndim from parts is the same as the user-defined NDIM above
        if(not dimData.NDIM == param.dict_eltype_ndim.get(elem_group.eltype)):
            print("WARNING: User-defined no. of dimensions is not compatible\
                  with that in the input file for eltype: "+elem_group.eltype)
        # check if nst from parts is the same as the user-defined NST above
        if(not dimData.NST == param.dict_eltype_nst.get(elem_group.eltype)):
            print("WARNING: User-defined no. of strains is not compatible\
                  with that in the input file for eltype: "+elem_group.eltype)
        # check if ndof_node from parts is the same as the user-defined NDOF_NODE above
        if(not dimData.NDOF_NODE == param.dict_eltype_ndofnode.get(elem_group.eltype)):
            print("WARNING: User-defined no. of dofs per node is not compatible\
                  with that in the input file for eltype: "+elem_group.eltype)
        # check if the element types in the input file are supported by the programme
        if(not elem_group.eltype in param.tuple_supported_eltypes):
            print("WARNING: eltype in the input file is not one of the \
                  supported element types! This element will not be used in \
                  the analysis: "+elem_group.eltype)
        # check if the element types in the input file are the intended ones of the user
        if(not elem_group.eltype in dimData.ELEM_TYPES):
            print("WARNING: eltype: in the input file is not one of the \
                  user-defined/expected element types! This element will not \
                  be used in the analysis: "+elem_group.eltype)
    # check ndof_node for elem_groups
    # NOTE: same no. of dofs per node for all elem groups in a part    
    def same_ndof_node(elem_groups):
        x = []
        for elem_group in elem_groups:
            x.append(param.dict_eltype_ndofnode.get(elem_group.eltype))
        return x.count(x[0]) == len(x)
    
    if(not same_ndof_node(part.elem_groups)):
        print("WARNING: Nodes in a part should have the same no. of dofs!")



def form_nodes(part):
    return np.asarray(part.nodes)



def form_elem_lists(part, NDOF_NODE, ELEM_TYPES, dict_elset_matID):
    """form the lists of elements for this part. They could be objects with 
    components and associated methods. requires user-defined no. of dofs per node 
    and expected element types for this part"""
    
    # inner function to form the cnc_dof from cnc_node
    def form_elem_cnc_dof(elem_cnc_node, ndof_node):
        """ form cnc_dof (dof connectivity) of the element based on 
        its cnc_node (nodal connectivity) and ndof_node (no. of dofs per node) """
        elem_cnc_dof = []
        for jnd in elem_cnc_node:
            # Python convension starts from 0 when setting the range below
            elem_cnc_dof.extend(list(range(jnd*ndof_node, (jnd+1)*ndof_node)))
        return elem_cnc_dof
    
    # inner function to update the matID of an element in the elem_lists object
    def update_elem_matID(elem_lists, ie, matID):
        istart = 0
        for elem_list in elem_lists:
            if istart <= ie < istart+len(elem_list.elems):
                elem_list.elems[ie].matID = matID
                break
            istart += len(elem_list.elems)
    
    
    elem_lists = []
    for elem_group in part.elem_groups:
        if elem_group.eltype in param.tuple_supported_eltypes and \
           elem_group.eltype in ELEM_TYPES:
            # initialize the empty list of elems for this eltype
            elems = []
            # loop over all elems in the nodal connectivity matrix
            for elem_cnc_node in elem_group.cnc_node:
                # obtain the element's dof connectivity list
                elem_cnc_dof = form_elem_cnc_dof(elem_cnc_node, NDOF_NODE)
                #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # form this element based on its eltype and using its node and 
                # dof cnc lists; to be updated with every new element class
                #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if elem_group.eltype in param.tuple_frame2d2_eltypes:
                    elems.append(frame2d2elem.frame2d2elem(elem_cnc_node, elem_cnc_dof))
                elif elem_group.eltype in param.tuple_tri2d3_eltypes:
                    elems.append(tri2d3elem.tri2d3elem(elem_cnc_node, elem_cnc_dof))
                else:
                    print('WARNING: this eltype is not added in elem_lists:'\
                          +elem_group.eltype)
            # form the elem_list for this eltype and append to the full elem_lists
            elem_lists.append(elem_list(elem_group.eltype, elems))
    
    # update matID of each element if dict_elset_matID it is not empty
    if dict_elset_matID:
        for elset in part.elsets:
            if dict_elset_matID.get(elset.name) is not None:
                matID = dict_elset_matID.get(elset.name)
                for ie in elset.setlist:
                    update_elem_matID(elem_lists, ie, matID)
    
    return elem_lists



def form_bcds_cloads(part, dict_nset_data, NDOF_NODE):
    """Form the boundary condition (bcd) and cload lists based on the data of part
    and the user-defined interpreter of nset names and user-expected no. of dofs 
    per node"""
    # initialize the lists
    bcd_dofs    = []
    bcd_values  = []
    cload_dofs   = []
    cload_values = []
    # form the lists
    for nset in part.nsets:
        nsetdata = dict_nset_data.get(nset.name)
        if nsetdata is None:
            print("WARNING: unrecognized nset name:"+nset.name)
        else:
            if nsetdata.settype == 'bcd':
                for inode in nset.setlist:
                    bcd_dofs += [inode*NDOF_NODE+x for x in nsetdata.nodalDofs]
                    bcd_values += nsetdata.dofValues
            elif nsetdata.settype == 'cload':
                for inode in nset.setlist:
                    cload_dofs += [inode*NDOF_NODE+x for x in nsetdata.nodalDofs]
                    cload_values += nsetdata.dofValues
            else:
                print("WARNING: unrecognized nsetdata type:"+nsetdata.settype)
        
    return [bcd_dofs, bcd_values, cload_dofs, cload_values]



def form_list_dload_data(part, dict_elset_dloadID, dload_functions):
    list_dload_data = []
    
    if dict_elset_dloadID:
        for elset in part.elsets:
            if dict_elset_dloadID.get(elset.name) is not None:
                dloadID = dict_elset_dloadID.get(elset.name)
                list_dload_data.append(dload_data(elset.setlist, dload_functions[dloadID]))
    
    return list_dload_data