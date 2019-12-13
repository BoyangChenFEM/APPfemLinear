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
jobname = 'Tri3fFEMdload'
inputfile = jobname+'.inp'

# Define dimensional parameters
NDIM = 2 # no. of dimensions
NST  = 3 # no. of strains/stresses
NDOF_NODE = 2 # no. of dofs per node
ELEM_TYPES = ('CPS3', 'CPE3') # expected element types in the inp file
dimData = dimension_data(NDIM, NST, NDOF_NODE, ELEM_TYPES)

# Define material parameters
E  = 200000.  # Young's modulus, MPa
nu = 0.3     # Poisson ratio
t  = 50.      # thickness out of plane, mm
is_planestress = True
# define the list of material sections
Materials = [linear_elastic.isotropic2D(E, nu, t, is_planestress)]
#------------------------------------------------------------------------------
# define the interpreter connecting elset names to Materials list index 
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dict_elset_matID = {}

# Define the concentrated loads and bcds
P = -1000. # N
applied_disp_u1 = 1 # mm
applied_disp_u2 = -3 # mm
#------------------------------------------------------------------------------
# Define interpretations of nset names to form lists of bcds and loads
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dict_nset_data = {\
'fix-all'   : nset_data(settype='bcd', nodalDofs=list(range(NDOF_NODE)), dofValues=[0]*NDOF_NODE),\
'fix-dir1'  : nset_data(settype='bcd', nodalDofs=[0], dofValues=[0]),\
'fix-dir2'  : nset_data(settype='bcd', nodalDofs=[1], dofValues=[0]),\
'disp-dir1' : nset_data(settype='bcd', nodalDofs=[0], dofValues=[applied_disp_u1]),\
'disp-dir2' : nset_data(settype='bcd', nodalDofs=[1], dofValues=[applied_disp_u2]),\
'cload-dir1': nset_data(settype='cload', nodalDofs=[0], dofValues=[P]),\
'cload-dir2': nset_data(settype='cload', nodalDofs=[1], dofValues=[P]) }

# Define the distributed loads (user-defined)
def func_dload(pos_vec):
    x = pos_vec[0]
    y = pos_vec[1]
    # define x component of the applied traction function
    tx = 0
    # define y component of the applied traction function
    eps = 1.E-9 # a small number for tolerance
    y_surf = 50
    if abs(y-y_surf) < eps:
        ty = -2 # N/mm
    return [tx, ty]
dload_functions = [dload_function(expression=func_dload, order=0)]

#------------------------------------------------------------------------------
# define the interpreter connecting elset names to dloads list index 
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dict_elset_dloadID = {'dload1':0}


###############################################################################
# Call the kernel programme (same for different jobs)
# return all data for postprocessing
###############################################################################
[InputPartsData, nodes, elem_lists, f, a, RF] = \
kernel_program(inputfile, dimData, Materials, dict_elset_matID, \
               dict_nset_data, dload_functions, dict_elset_dloadID)




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
vtkoutput_ig(jobname, elem_lists, NST)
