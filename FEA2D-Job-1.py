"""
Created on Fri Nov 01 2019
Finite Element Programme
Job file for linear elastic analysis under concentrated loads
The main program reads mesh and nsets from the part definition of Abaqus input 
file (only a single part is supported) and form the lists of nodes and elems.
It forms the lists of bcds and cloads (concentrated loads) based on user-
defined nset name reader 'READ_NSET_NAME' 
@author: Boyang CHEN TU Delft
"""

# import the main program interface of LinearFEMProgram
from LinearFEMProgram import Kernel
# import the default output program to write outputs in vtk format
from LinearFEMProgram.output import vtkoutput1part as vtkoutput, \
                                    vtkoutput1part_igpoints as vtkoutput_ig
## for debugging
#from LinearFEMProgram.Preprocessing import read_abaqus_parts as read_abaqus
#from LinearFEMProgram.Kernel import form_nodes_elems, form_bcds_cloads
#import numpy as np
## end debugging

###############################################################################
# User-defined problem-specific parameters
# for dimension, no. of strains, no. of dofs per node, material parameters,
# material stiffness matrix, load values, and nset names of the applied 
# loads and boundary conditions
###############################################################################

# import the right material module
from LinearFEMProgram.Materials import linear_elastic
# import modules for user-defined fast visualization in Python console
import matplotlib.pyplot as plt


# set input file name
jobname = 'salome-test3disp'
#jobname = 'Tri3fFEMmesh1'
inputfile = jobname+'.inp'

# Define dimensional parameters
NDIM = 2 # no. of dimensions
NST  = 3 # no. of strains/stresses
NDOF_NODE = 2 # no. of dofs per node
ELEM_TYPES = ('CPS3', 'CPE3') # expected element types in the inp file

# Define material parameters
E  = 200000.  # Young's modulus, MPa
nu = 0.3     # Poisson ratio
t  = 50.      # thickness out of plane, mm
is_planestress = True
# define the material stiffness matrix
Dmat = t*linear_elastic.Dmatrix_2D_isotropic(E, nu, is_planestress)

# Define the load
P = -1000. # N
applied_disp_u1 = 1 # mm
applied_disp_u2 = -3 # mm

#------------------------------------------------------------------------------
# Define interpretations of nset names to form lists of bcds and loads
# This will be called in the kernel programme
#------------------------------------------------------------------------------
def READ_NSET_NAME(set_name, node_list, ndof_node):
    """A function to get the lists of bcd/load dofs and values from nset name,
    node list and ndof_node""" 
    bcd_dofs    = []
    bcd_values  = []
    cload_dofs   = []
    cload_values = []
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if set_name.lower() == 'fix-all': # fully clamped/encastre
        for inode in node_list:
            bcd_dofs   += range(inode*ndof_node, inode*ndof_node+ndof_node)
            bcd_values += [0]*ndof_node
    elif set_name.lower() == 'fix-dir1':
        for inode in node_list:
            bcd_dofs   += [inode*ndof_node]
            bcd_values += [0]
    elif set_name.lower() == 'fix-dir2':
        for inode in node_list:
            bcd_dofs   += [inode*ndof_node+1]
            bcd_values += [0]
    elif set_name.lower() == 'disp-dir1':
        for inode in node_list:
            bcd_dofs   += [inode*ndof_node]
            bcd_values += [applied_disp_u1]
    elif set_name.lower() == 'disp-dir2':
        for inode in node_list:
            bcd_dofs   += [inode*ndof_node+1]
            bcd_values += [applied_disp_u2]
    elif set_name.lower() == 'cload-dir1': # concentrated load along direction 1
        for inode in node_list:
            cload_dofs   += [inode*ndof_node]
            cload_values += [P]
    elif set_name.lower() == 'cload-dir2': # concentrated load along direction 2
        for inode in node_list:
            cload_dofs   += [inode*ndof_node+1]
            cload_values += [P]
    else:
        #raise ValueError("BCD/cload nset name not supported!")
        print("WARNING: BCD/cload nset name "+set_name+" is not supported!")
        
    return [bcd_dofs, bcd_values, cload_dofs, cload_values]


## debugging 
#parts = read_abaqus.read_parts_from_inputfile(inputfile)
##print(parts[0].name)
##print(parts[0].elem_groups[0].eltype)
#print(parts[0].elem_groups[0].cnc_node)
#print(parts[0].nsets[1].setlist)
#print(parts[0].nsets[2].setlist)
#nodes, elems = form_nodes_elems(parts[0])
#[bcd_dofs, bcd_values, cload_dofs, cload_values] = \
#form_bcds_cloads(parts[0], READ_NSET_NAME)
#print(elems[0].cnc_node)
#print(elems[0].cnc_dof)
##print(nodes[elems[0].cnc_node])
#print(bcd_dofs)
#print(bcd_values)
#print(cload_dofs)
#print(cload_values)
## end debugging


###############################################################################
# Call the kernel programme (same for different jobs)
# return all data for postprocessing
###############################################################################
[InputPartsData, nodes, elem_lists, f, a, RF] = \
Kernel.kernel_program(inputfile, NDIM, NST, NDOF_NODE, ELEM_TYPES, Dmat, \
                      READ_NSET_NAME)




###############################################################################
# Postprocessing (Problem-specific) 
###############################################################################
## Primary output
##print('Nodal displacement at bcd points (mm):')
##print(a[bcd_dofs])
#print('Nodal displacement at load point (mm):')
#print(a[cload_dofs])
#print('RF at bcd points (N):')
#print(RF[bcd_dofs])
##print('RF at load point (mm):')
##print(RF[cload_dofs])
#

#------------------------------------------------------------------------------
# simple plotting within Python console for fast visualization
#------------------------------------------------------------------------------
r=100 # scaling factor for deformed plot
#plt.figure(figsize=[6.4,9.6])
# plot undeformed nodes
plt.plot(nodes[:,0],nodes[:,1],'.b',label='initial mesh')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.legend(loc='best')
#plt.savefig("initialmesh.png", dpi=250)
# plot deformed nodes
x_curr = nodes[:,0] + r*a[0::NDOF_NODE].reshape(len(nodes))
y_curr = nodes[:,1] + r*a[1::NDOF_NODE].reshape(len(nodes))
#plt.scatter(x_curr,y_curr,c='r',alpha=0.5,label='deformed')
plt.plot(x_curr, y_curr, '.r', label='deformed mesh (scale=%d)' %r)
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.legend(loc='best')
#plt.savefig("deformedmeshscale100.pdf")
plt.show()


#------------------------------------------------------------------------------
# Default VTK output of data for visualization using Paraview
#------------------------------------------------------------------------------
vtkoutput(jobname, nodes, elem_lists, f, a, RF, NDIM, NDOF_NODE)
vtkoutput_ig(jobname, elem_lists, NST)
          

#
#print('Reaction forces (N) are:')
#print(RF)

## Derive strain & stress for frame elements
#strain = np.zeros([nelems,NST])
#stress = np.zeros([nelems,NST])
#for e in range(nelems):
#    # extract the nodal coordinates of the element
#    elem_nodes = node_coords[elem_cnc[e,:],:]
#    # extract the nodal printlacements of the element
#    elem_dofs  = a[elem_cnc_dof[e,:]]
#    # element-specific calculations of strain and stress below
#    [strain[e,:],stress[e,:]] = tri2d3n.strain_stress(elem_nodes, Dmat, elem_dofs)