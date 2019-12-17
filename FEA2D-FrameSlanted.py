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
# import the element stress-strain function for postprocessing
from LinearFEMProgram.Elements import frame2d2elem
# import modules for user-defined fast visualization in Python console
import matplotlib.pyplot as plt
import numpy as np
from prettytable import PrettyTable


# set input file name
jobname = 'SlantedFrame1'
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
A_frame = b*h    # mm^2
I       = b*h**3/12 # mm^4

# define the list of materials/sections
Materials = [linear_elastic.frame2D(E, A_frame, I)]
#------------------------------------------------------------------------------
# define the interpreter connecting elset names to Materials list index 
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# dict_elset_matID = {}

#------------------------------------------------------------------------------
# Define interpretations of nset names to form lists of bcds and loads
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dict_nset_data = {\
'Left'   : nset_data(settype='bcd', nodalDofs=[0, 1, 2], dofValues=[0, 0, 0]),\
'Right'  : nset_data(settype='bcd', nodalDofs=[0, 1, 2], dofValues=[0, 0, 0]) }

# Define the distributed loads (user-defined)
# for this problem, it is defined in frame's local coordinates, linearly-varied,
# orineted along the -Y(global) direction.
def func_linear_loc(xloc):
    # define geometrical parameters to define the local frame orientation 
    L = 2000.0
    H = 1500.0
    S = 2500.0
    costheta = L/S
    sintheta = H/S
    # ref. point local coords, where the curve parameter s = 0 (start of line)
    refP = [0, 0]
    # max. amplitude of the distribution along the frame (force per frame length)
    q0 = 1.6 #N/mm    
    # define local x component of the applied traction function
    tx = -q0*(xloc[0]-refP[0])/S*sintheta
    # define local y component of the applied traction function
    ty = -q0*(xloc[0]-refP[0])/S*costheta
    return [tx, ty]

# the following function defined on global x-y coords (for both inputs and outputs)
# produces the same result as the local one above. The key aspect is that this 
# function must be defined on the boundary itself, not along any projection of it;
# if it is the latter, then it must be transformed to a distribution along the 
# boundary while taking into account of the change of domain, i.e., the jacobian.
def func_linear_glb(xg):
    # define geometrical parameters to define the local frame orientation 
    L = 2000.0
    H = 1500.0
    # ref. point global coords, where the curve parameter s = 0 (start of line)
    refP = [0, 0]
    # max. amplitude of the distribution along the frame (force per frame length)
    q0 = 1.6 #N/mm    
    # tolerance for zero, 0.1 mm
    tol = 1.e-1 
    # initialize tx and ty
    tx = 0
    ty = 0
    # apply traction on the line
    if abs((xg[1]-refP[1])*L - (xg[0]-refP[0])*H) < tol:
        # define x component of the applied traction function
        tx = 0
        # define y component of the applied traction function
        ty = -q0*xg[0]/L
    return [tx, ty]

# 'local' indicates that tx and ty are along local frame axis and dload expression
# takes in point coordinates in local coordinate system
# order=0 indicates that the loading is constant
# order=1 indicates that the loading is linearly-varied
#dload_functions = [dload_function(expression=func_linear_loc, coord_system='local', order=1)]
dload_functions = [dload_function(expression=func_linear_glb, coord_system='global', order=1)]
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
               dict_elset_dload=dict_elset_dload)




###############################################################################
# Postprocessing (Problem-specific) 
###############################################################################
#------------------------------------------------------------------------------
# Default VTK output of data for visualization using Paraview
#------------------------------------------------------------------------------
vtkoutput(jobname, nodes, elem_lists, f, a, RF, NDIM, NDOF_NODE)

# simple printing of f and RF
print(f)
print(RF)

#------------------------------------------------------------------------------
# simple plotting within Python console for visualization of mesh
#------------------------------------------------------------------------------
r=100 # scaling factor for deformed plot
#plt.figure(figsize=[6.4,9.6])
plt.figure()
for elist in elem_lists:
    for elem in elist.elems:
        elnodes = nodes[elem.cnc_node]
        eldofs  = a[elem.cnc_dof]
        # undeformed
        x_init = elnodes[:,0]
        y_init = elnodes[:,1]
        p1, =plt.plot(x_init,y_init,'.b-')
        # deformed
        x_curr = elnodes[:,0] + r*eldofs[0::NDOF_NODE].reshape(len(elnodes))
        y_curr = elnodes[:,1] + r*eldofs[1::NDOF_NODE].reshape(len(elnodes))
        p2, =plt.plot(x_curr, y_curr, '.r-')
p1.set_label('undeformed')
p2.set_label('deformed(scale=%d)' %r)
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.legend(loc='best')
plt.show()
#plt.savefig("deformedmeshscale100.png", dpi=250)
#plt.savefig("deformedmeshscale100.pdf")

# Derive strain & stress for frame elements
allstrain=[]
allstress=[]
for elist in elem_lists:
    for elem in elist.elems:
        elnodes = nodes[elem.cnc_node]
        eldofs  = a[elem.cnc_dof]
        # element-specific calculations of strain and stress below
        # max stress occurs at 1 of the 4 corners
        # left end relative coordinate: 0
        # right end relative coordinate: 1
        # top surf relative coordinate: 1/2
        # bot surf relative coordinate: -1/2
        points = [[0,-1/2],[0,1/2],[1,-1/2],[1,1/2]]
        elstress = []
        elstrain = [] 
        for point in points:
            eps, sig = frame2d2elem.strain_stress(elnodes, eldofs, \
                                            Materials[elem.matID], h, point)
            elstrain.append(eps.tolist())
            elstress.append(sig.tolist())
        allstrain.append(elstrain)
        allstress.append(elstress)
        


# Make tablular data for nodes
af = a.reshape(len(nodes),NDOF_NODE)
af = np.around(af,decimals=4)

tab1 = PrettyTable()
tab1.field_names=["Nodes", "u (mm)", "v (mm)", "theta (rad) for frame nodes"]
for n in range(len(nodes)):
    tab1.add_row([n+1, af[n,0], af[n,1], af[n,2]])

print(tab1)
#with open('tab1.txt','w') as tab1txt:
#    tab1txt.write(str(tab1))

# Make tabular data for elements
tab2 = PrettyTable()
tab2.field_names=["Elements", "Max. absolute axial stress (MPa)"]
for e, elstress in enumerate(allstress):
    elstress = sum(elstress, [])
    maxabsstress = abs(max(elstress, key=abs))
    tab2.add_row([e+1, round(maxabsstress,3)])

print(tab2)
#with open('tab2.txt','w') as tab2txt:
#    tab2txt.write(str(tab2))
