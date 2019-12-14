"""
Created on Fri Nov 01 2019
Finite Element Programme
Job file for linear elastic analysis under concentrated loads
The main program reads mesh and nsets from the part definition of Abaqus input 
file (only a single part is supported) and form the lists of nodes and elems.
It forms the lists of bcds, cloads (concentrated loads) and dloads (distributed
loads). It supports multiple materials and multiple element types with the 
same number of nodal DoFs: 
- tri3, quad4, tri6 and quad8 are compatible
- frame2d and truss2d are not compatible
- truss must be modelled by frame elems with a truss stiff matrix
@author: Boyang CHEN TU Delft
"""

# import the main program interface of LinearFEMProgram
from LinearFEMProgram.Kernel import dimension_data, nset_data, dload_function,\
                                    kernel_program
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
#jobname = 'salome-test3disp'
jobname = 'FrameTruss'
inputfile = jobname+'.inp'

# Define dimensional parameters
NDIM = 2 # no. of dimensions
NST  = 1 # no. of strains/stresses
NDOF_NODE = 3 # no. of dofs per node
ELEM_TYPES = ('T2D2', 'B23') # expected element types in the inp file
dimData = dimension_data(NDIM, NST, NDOF_NODE, ELEM_TYPES)

# Define material & section parameters
E       = 10000.  # Young's modulus, MPa
b       = 100    # mm
h       = 100    # mm
A_truss = 1000   # mm^2
A_frame = b*h    # mm^2
I       = b*h**3/12 # mm^4

# define the list of material sections
Materials = [linear_elastic.truss(E, A_truss), linear_elastic.frame2D(E, A_frame, I)]
#------------------------------------------------------------------------------
# define the interpreter connecting elset names to Materials list index 
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dict_elset_matID = {'TrussElems':0, 'FrameElems':1}

# Define the concentrated loads and bcds
P = -5000. # N
#------------------------------------------------------------------------------
# Define interpretations of nset names to form lists of bcds and loads
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dict_nset_data = {\
'FrameBL'   : nset_data(settype='bcd', nodalDofs=list(range(NDOF_NODE)), dofValues=[0]*NDOF_NODE),\
'FrameTR'   : nset_data(settype='bcd', nodalDofs=list(range(NDOF_NODE)), dofValues=[0]*NDOF_NODE),\
'TrussBR'   : nset_data(settype='bcd', nodalDofs=[0, 1], dofValues=[0, 0]),\
'LoadPoint' : nset_data(settype='cload', nodalDofs=[1], dofValues=[P]) }

# Define the distributed loads (user-defined)
# for this problem, it is defined in frame's local coordinates, uniformly
# distributed.
def func_dload():
    # define local x component of the applied traction function
    tx = 0
    # define local y component of the applied traction function
    ty = -1 # N/mm
    return [tx, ty]
dload_functions = [dload_function(expression=func_dload, coord_system='local', order=0)]

#------------------------------------------------------------------------------
# define the interpreter connecting elset names to dload functions above 
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dict_elset_dload = {'q_elems': dload_functions[0]}


###############################################################################
# Call the kernel programme (same for different jobs)
# return all data for postprocessing
###############################################################################
[InputPartsData, nodes, elem_lists, f, a, RF] = \
kernel_program(inputfile, dimData, Materials, dict_elset_matID, \
               dict_nset_data, dict_elset_dload)




###############################################################################
# Postprocessing (Problem-specific) 
###############################################################################
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
