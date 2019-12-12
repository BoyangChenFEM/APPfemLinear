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
from LinearFEMProgram.Kernel import nset_data, kernel_program
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
jobname = 'Tri3fFEMmesh1'
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




###############################################################################
# Call the kernel programme (same for different jobs)
# return all data for postprocessing
###############################################################################
[InputPartsData, nodes, elem_lists, f, a, RF] = \
kernel_program(inputfile, NDIM, NST, NDOF_NODE, ELEM_TYPES, Dmat, dict_nset_data)




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
