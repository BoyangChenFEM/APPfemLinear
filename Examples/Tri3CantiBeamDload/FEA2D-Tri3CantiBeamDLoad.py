"""
Job script for linear elastic analysis

for a rectangular cantilever beam under end force loading, meshed with linear
triangular elements.

@author: Boyang CHEN, TU Delft Aerospace, 2019
"""

# import the main program interface of LinearFEMProgram
from LinearFEMProgram.Kernel import dimension_data, nset_data, dload_function,\
                                    kernel_program
# import the default output program to write outputs in vtk format
from LinearFEMProgram.output import vtkoutput1part as vtkoutput

###############################################################################
# User-defined problem-specific parameters
###############################################################################

# import the right material module
from LinearFEMProgram.Materials import linear_elastic
# import program to write igpoint outputs in vtk format (for solid elems only)
from LinearFEMProgram.output import vtkoutput1part_igpoints as vtkoutput_ig
# import modules for user-defined fast visualization in Python console
import matplotlib.pyplot as plt
import numpy as np


# set input file name
jobname = 'Tri3fFEMdload'
inputfile = jobname+'.inp'

# Define dimensional parameters
NDIM = 2 # no. of dimensions
NDOF_NODE = 2 # no. of dofs per node
ELEM_TYPES = ('CPS3', 'CPE3') # expected element types in the inp file
dimData = dimension_data(NDIM, NDOF_NODE, ELEM_TYPES)

# Define material parameters
E  = 200000.  # Young's modulus, MPa
nu = 0.3     # Poisson ratio
t  = 50.      # thickness out of plane, mm
is_planestress = True
# define the list of materials/sections
Materials = [linear_elastic.isotropic2D(E, nu, t, is_planestress)]

# Define the concentrated loads and bcds

#------------------------------------------------------------------------------
# Define interpretations of nset names to form lists of bcds and loads
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# nodes in nset named 'fix-all' will be fully constrained
# nodes in nset named 'cload-dir2' will be loaded in force along y by
# a value represented by the variable P
dict_nset_data = {\
'fix-all'   : nset_data(settype='bcd', nodalDofs=list(range(NDOF_NODE)), dofValues=[0]*NDOF_NODE)}

# Define the distributed loads (user-defined)
# for this problem, it is defined in frame's local coordinates, uniformly
# distributed.
def func_dload(coords):
    # define local x component of the applied traction function
    tx = 0 # no axial loading
    # define local y component of the applied traction function
    # eps = 1.e-9 # small number as tolerance for zero
    # define applied vertical traction on top surf of beam (y=50)
    #if abs(coords[1]-50.0)<eps:
    if coords[1] == 50.0:
        ty = -1*t # N/mm # uniform vertical loading
    else:
        ty = 0
    return [tx, ty]

def func_dload1(coords):
    # define local x component of the applied traction function
    tx = 0 # no axial loading
    # define local y component of the applied traction function
    # eps = 1.e-9 # small number as tolerance for zero
    # define applied vertical traction on top surf of beam (y=50)
    #if abs(coords[1]-50.0)<eps:
    if coords[1] == 50.0:
        ty = (coords[0]-350.0)/1000.0*t # N/mm # linearly-varying vertical loading
    else:
        ty = 0
    return [tx, ty]
# 'local' indicates that tx and ty are along local frame axis 
# order=0 indicates that the loading is constant
dload_functions = [dload_function(expression=func_dload, coord_system='global', order=0),\
                   dload_function(expression=func_dload1, coord_system='global', order=1)]

#------------------------------------------------------------------------------
# define the interpreter connecting elset names to dload functions above 
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# all elements in the elset named 'q_elems' will be subjected to the distributed
# loading defined in dload_functions[0]
dict_elset_dload = {'dload1': dload_functions[1]}

###############################################################################
# Call the kernel programme (same for different jobs)
# return all data for postprocessing
###############################################################################
[InputPartsData, nodes, elem_lists, f, a, RF] = \
kernel_program(inputfile, dimData, Materials, dict_nset_data, \
               dict_elset_dload=dict_elset_dload)




###############################################################################
# Postprocessing (Problem-specific) 
###############################################################################
#------------------------------------------------------------------------------
# Default VTK output of data for visualization using Paraview
#------------------------------------------------------------------------------
vtkoutput(jobname, nodes, elem_lists, f, a, RF, NDIM, NDOF_NODE)

##------------------------------------------------------------------------------
## simple plotting within Python console for visualization of mesh
##------------------------------------------------------------------------------
#r=10 # scaling factor for deformed plot
##plt.figure(figsize=[6.4,9.6])
#plt.figure()
#for elist in elem_lists:
#    for elem in elist.elems:
#        elnodes = nodes[elem.cnc_node]
#        eldofs  = a[elem.cnc_dof]
#        # undeformed
#        x_init = np.append(elnodes,[elnodes[0]], axis=0)[:,0]
#        y_init = np.append(elnodes,[elnodes[0]], axis=0)[:,1]
#        p1, =plt.plot(x_init,y_init,'.b-')
#        # deformed
#        x_curr = elnodes[:,0] + r*eldofs[0::NDOF_NODE].reshape(len(elnodes))
#        y_curr = elnodes[:,1] + r*eldofs[1::NDOF_NODE].reshape(len(elnodes))
#        x_curr = np.append(x_curr, x_curr[0])
#        y_curr = np.append(y_curr, y_curr[0])
#        p2, =plt.plot(x_curr, y_curr, '.r-')
#p1.set_label('undeformed')
#p2.set_label('deformed(scale=%d)' %r)
#plt.xlabel('x (mm)')
#plt.ylabel('y (mm)')
#plt.legend(loc='best')
#plt.show()
##plt.savefig("deformedmeshscale100.png", dpi=250)
##plt.savefig("deformedmeshscale100.pdf")


#------------------------------------------------------------------------------
# Update element igpoints for output/postprocessing
# modify the x, u, strain and stress of integration points of each element
#------------------------------------------------------------------------------ 
for elist in elem_lists:
    for elem in elist.elems:
        elnodes = nodes[elem.cnc_node]
        eldofs  = a[elem.cnc_dof]
        elem.update_igpoints(elnodes, Materials[elem.matID], eldofs)
vtkoutput_ig(jobname, elem_lists)
