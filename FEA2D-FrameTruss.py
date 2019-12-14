"""
Job script for linear elastic analysis

for a combined frame-truss structure under combined distributed and concentrated
loadings. Hence, it has multiple material sections. Truss is modelled by frame
element with a truss material section.

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
# import modules for user-defined fast visualization in Python console
import matplotlib.pyplot as plt


# set input file name
jobname = 'FrameTruss'
inputfile = jobname+'.inp'

# Define dimensional parameters
NDIM = 2 # no. of dimensions
NDOF_NODE = 3 # no. of dofs per node
ELEM_TYPES = ('B23') # expected element types in the inp file
dimData = dimension_data(NDIM, NDOF_NODE, ELEM_TYPES)

# Define material & section parameters
E       = 10000.  # Young's modulus, MPa
b       = 100    # mm
h       = 100    # mm
A_truss = 1000   # mm^2
A_frame = b*h    # mm^2
I       = b*h**3/12 # mm^4

# define the list of materials/sections
Materials = [linear_elastic.truss(E, A_truss), linear_elastic.frame2D(E, A_frame, I)]
#------------------------------------------------------------------------------
# define the interpreter connecting elset names to Materials list index 
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Elements in the elset named 'TrussElems' will be attributed with Materials[0]
# Elements in the elset named 'FrameElems' will be attributed with Materials[1]
dict_elset_matID = {'TrussElems':0, 'FrameElems':1}

# Define the concentrated loads and bcds
P = -5000. # N
#------------------------------------------------------------------------------
# Define interpretations of nset names to form lists of bcds and loads
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Nset FrameBL & FrameTR: fully-constrained
# Nset TrussbR: pinned
# Nset LoadPoint: loaded with concentrated force along Y with value P
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
    tx = 0 # no axial loading
    # define local y component of the applied traction function
    ty = -1 # N/mm # uniform vertical loading
    return [tx, ty]
# 'local' indicates that tx and ty are along local frame axis 
# order=0 indicates that the loading is constant
dload_functions = [dload_function(expression=func_dload, coord_system='local', order=0)]

#------------------------------------------------------------------------------
# define the interpreter connecting elset names to dload functions above 
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# all elements in the elset named 'q_elems' will be subjected to the distributed
# loading defined in dload_functions[0]
dict_elset_dload = {'q_elems': dload_functions[0]}


###############################################################################
# Call the kernel programme (same for different jobs)
# return all data for postprocessing
###############################################################################
[InputPartsData, nodes, elem_lists, f, a, RF] = \
kernel_program(inputfile, dimData, Materials, dict_nset_data, \
               dict_elset_matID, dict_elset_dload)




###############################################################################
# Postprocessing (Problem-specific) 
###############################################################################
#------------------------------------------------------------------------------
# Default VTK output of data for visualization using Paraview
#------------------------------------------------------------------------------
vtkoutput(jobname, nodes, elem_lists, f, a, RF, NDIM, NDOF_NODE)

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


