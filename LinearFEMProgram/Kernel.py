"""
Core Program for linear elastic Finite Element Analysis

This program reads only the mesh, nsets and elsets from the part definition of 
Abaqus input file (only a single part is supported at the moment). 

All information of bcds, loads and materials (in case there are mutliple) must
be defined in the job script, connected to the nodes and elems in the part through
the nsets and elsets in the inputfile. The main program then forms the datalists 
of nodes, elems, bcds, cloads (concentrated loads) and dloads (distributed loads).

The program supports multiple materials. 

The program supports multiple element types, as long as they have the 
same number of nodal DoFs and their meanings are consistent. For example: 
- tri3, quad4, tri6 and quad8 are compatible
- frame2d and truss2d are not compatible; truss must be modelled by frame 
elements with a truss material section (supported by the frame element)

@author: Boyang CHEN, TU Delft Aerospace, 2019
"""
import numpy as np
from . import Parameters as param
from .Preprocessing import read_abaqus_parts as read_abaqus
from .Elements import frame2d2elem, tri2d3elem
from .Solvers import assembler, direct_solver_edu as solver

###############################################################################
###############################################################################
# Define classes for objects in kernel program
###############################################################################
###############################################################################
class dimension_data:
    """user-input dimensional data in the job script, used to verify with the 
    input file data. Its components are:
        - NDIM: number of dimensions in analysis space
        - NDOF_NODE: number of DoFs per node
        - ELEM_TYPES: expected element types in the model"""
    def __init__(self, NDIM, NDOF_NODE, ELEM_TYPES):
        self.NDIM = NDIM
        self.NDOF_NODE = NDOF_NODE
        self.ELEM_TYPES = ELEM_TYPES

class nset_data:
    """a datatype defined to store nodeset data for boundary conditions and 
    concentrated loads applied on the nodes. It has the following components:
        - settype: either 'bcd' for imposed displacement or 'cload' for applied
        concentrated load
        - nodalDofs: the Dof indices affected by the imposed bcd or cload
        - dofValues: the value of the applied displacement or cload"""    
    def __init__(self, settype, nodalDofs, dofValues):
        self.settype = settype # either 'bcd' or 'cload'
        self.nodalDofs = nodalDofs # indices of affected nodal DoFs
        self.dofValues = dofValues # values of imposed disp or cload
        
class dload_function:
    """a datatype defined to represent a generic distributed loading. Its 
    components are:
        - expression: the loading, a function of spacial coordinates
        - coord_system: coordinate system for the loading function's inputs & 
        outputs. If it is 'local', then the loading is defined in local
        coordinate system, as in the case of beam. If it is 'global', then the 
        loading is defined in the global coordinate system, as in solid elements.
        default: global
        - order: the order of the function, based on this the number of gauss 
        points will be determined when integrating the force vector. default=0,
        i.e., a constant distribution.
        - """
    def __init__(self, expression, coord_system='global', order=0):
        self.expression = expression
        self.coord_system = coord_system
        self.order = order

class dload_data:
    """ A dataset which contains:
        - elset: the element set that is applied with a distributed load
        - function: the dload_function object which represents the distributed 
        load that is applied on the elements in elset"""
    def __init__(self, elset, function):
        self.elset = elset
        self.function = function

class elem_list:
    """A list of elements of the same type.
    Its components are: 
    - eltype: the element type for this list of elements
    - elems: the indices of all the elements in the list"""
    def __init__(self, eltype, elems=None):
        self.eltype = eltype
        self.elems  = elems



###############################################################################
###############################################################################
#---------------------- The kernel program -----------------------------------#
###############################################################################
###############################################################################
def kernel_program(inputfile, dimData, Materials, dict_nset_data, \
                   dict_elset_matID={}, dict_elset_dload={}):
    """The kernel_program should be called by the job script (e.g., Job-1.py) 
    where the user defines: 
    - inputfile: the name of the input file
    - dimData: the dimensional data (see class dimension_data)
    - Materials: the list of materials used in the analysis (see package Elements)
    - dict_nset_data: the dictionary of nset_data (for bcds and concentrated loads)
    where the keys are nset names read from inputfile and values are nset_data as
    defined in the class nset_data
    and optionally: 
    - dict_elset_matID: a dictionary where each key is an elset name defined in 
    inputfile, and its value is the corresponding index of material in the Materials
    list for elements in this elset. This dictionary needs to be defined when 
    multiple materials/material sections are present in the model
    - dict_elset_dload: a dictionary where each key is an elset name defined in 
    inputfile, and its value is the corresponding dload_data (see class dload_data)
    for all elements in this elset, meaning that these elements are subjected to 
    the distributed loading defined by this dload_data. This is needed when 
    distributed loading is present in the model"""
    ###########################################################################
    # Preprocessing 
    ###########################################################################
    # Read data from Abaqus input file and form abaqus parts
    parts = read_abaqus.read_parts_from_inputfile(inputfile)
    
    # check if there is only one part
    # in the future, consider making a loop over all parts
    if(not len(parts)==1):
        raise ValueError('Only a single part is supported!')
    
    # verification of dimensional parameters before proceeding
    verify_dimensional_parameters(parts[0], dimData)
        
    # form lists of nodes and elem_lists (eltype and elem indices of this type)
    nodes = form_nodes(parts[0])
    elem_lists = form_elem_lists(parts[0], dimData.NDOF_NODE, dimData.ELEM_TYPES,\
                                 dict_elset_matID)
    
    # form lists of bcds and cloads
    [bcd_dofs, bcd_values, cload_dofs, cload_values] = \
    form_bcds_cloads(parts[0], dict_nset_data, dimData.NDOF_NODE)
    
    # form lists of elset for distributed loads
    list_dload_data = form_list_dload_data(parts[0], dict_elset_dload)

    
    ###########################################################################
    # Assembler 
    # obtain the full stiffness matrix K and external distributed force vector f
    ###########################################################################
    # form the list of all the elems for assembly
    elems = []
    for elist in elem_lists:
        elems.extend(elist.elems)
    # call assembler
    [K, f] = assembler(nodes, elems, dimData.NDOF_NODE, Materials, list_dload_data)
    
    
    ###########################################################################
    # Solver
    # modify the stiffness matrix and force vector based on applied bcds and loads
    # obtain dof vector a and reaction force vector RF, both size ndof by 1
    ###########################################################################   
    [a, RF] = solver(K, f, bcd_dofs, bcd_values, cload_dofs, cload_values)
    
    
    return [parts, nodes, elem_lists, f, a, RF]







###############################################################################
###############################################################################
# Auxilliary programs
###############################################################################
###############################################################################
    
def verify_dimensional_parameters(part, dimData):
    """A function to perform sanity checks at the meta level before analysis"""
    
    for elem_group in part.elem_groups:
        # check if ndim from parts is the same as the user-defined NDIM above
        if(not dimData.NDIM == param.dict_eltype_ndim.get(elem_group.eltype)):
            print("WARNING: User-defined no. of dimensions is not compatible\
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
            print("WARNING: an eltype in the input file is not one of the \
                  user-defined/expected element types! This element will not \
                  be used in the analysis: "+elem_group.eltype)
    
    def same_ndof_node(elem_groups):
        """an inner function to check ndof_node for elem_groups; there should be 
        the same no. of dofs per node for all elem groups in a part"""
        x = []
        for elem_group in elem_groups:
            x.append(param.dict_eltype_ndofnode.get(elem_group.eltype))
        return x.count(x[0]) == len(x)
    
    if(not same_ndof_node(part.elem_groups)):
        print("WARNING: Nodes in a part should have the same no. of dofs!")



def form_nodes(part):
    """change the list into numpy arrays for better numerical operations"""
    return np.asarray(part.nodes)



def form_elem_lists(part, NDOF_NODE, ELEM_TYPES, dict_elset_matID={}):
    """form the lists of elements for this part. They could be objects with 
    components and associated methods. requires user-defined no. of dofs per node 
    and expected element types for this part; in case multiple materials are 
    present, also dict_elset_matID which relates elset with their material ID
    in the Materials list"""
    
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
                elem_list.elems[ie-istart].matID = matID
                break
            istart += len(elem_list.elems)
    
    # form the list of elem_list object (see elem_list class above) for this model
    # basically, it is a data structure which groups elements of the same type 
    # together first into an elem_list object, then group all the elem_list objects
    # into a list of element_list objects: elem_lists (not the best name, I know)
    elem_lists = []
    for elem_group in part.elem_groups:
        if elem_group.eltype in param.tuple_supported_eltypes and \
           elem_group.eltype in ELEM_TYPES:
            # initialize the empty list of elems for this eltype
            elems = []
            # loop over all rows of the nodal connectivity matrix
            # each row is basically the nodal connectivity of an element
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
    if dict_elset_matID:# true of it is not empty
        for elset in part.elsets:
            if dict_elset_matID.get(elset.name) is not None:
                # extract the material ID from the user-defined interpreter of
                # elset name passed in from the job script
                matID = dict_elset_matID.get(elset.name)
                # update the material ID to all the elements in this elset
                for ie in elset.setlist:
                    update_elem_matID(elem_lists, ie, matID)
    
    return elem_lists



def form_bcds_cloads(part, dict_nset_data, NDOF_NODE):
    """Form the boundary condition (bcd) and cload lists based on 
    - the data of part
    - the user-defined interpreter of nset names (dict_nset_data) 
    - no. of dofs per node (NDOF_NODE)"""
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



def form_list_dload_data(part, dict_elset_dload):
    """Read the elset name and elem indices in part, and use elset name to 
    extract dload function from dict_elset_dload provided by the user, and 
    form a list of dataset, each dataset is a pair like the following: 
    (elset, dload function), represented by the class 'dload_data'"""
    # initialize the list
    list_dload_data = []
    # form the list
    if dict_elset_dload:# true of it is not empty
        for elset in part.elsets:
            if dict_elset_dload.get(elset.name) is not None:
                # extract the dload function from the user-defined interpreter of
                # the elset name passed in from the job script
                dload_function = dict_elset_dload.get(elset.name)
                # group the elset elements and the dload function in the object
                # dload_data
                list_dload_data.append(dload_data(elset.setlist, dload_function))
    
    return list_dload_data