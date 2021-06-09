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
from LinearFEMProgram.Elements import truss2d2elem, frame2d2elem
# import modules for user-defined fast visualization in Python console
import LinearFEMProgram.Parameters as param
import matplotlib.pyplot as plt
import numpy as np
from prettytable import PrettyTable


# set input file name
jobname = 'FrameTruss'
inputfile = jobname+'.inp'

# Define dimensional parameters
NDIM = 2 # no. of dimensions
NDOF_NODE = 3 # max. no. of dofs per node
ELEM_TYPES = ('T2D2','B23') # expected element types in the inp file
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
'FrameBL'   : nset_data(settype='bcd', nodalDofs=[0,1,2], dofValues=[0,0,0]),\
'FrameTR'   : nset_data(settype='bcd', nodalDofs=[0,1,2], dofValues=[0,0,0]),\
'TrussBR'   : nset_data(settype='bcd', nodalDofs=[0,1], dofValues=[0,0]),\
'LoadPoint' : nset_data(settype='cload', nodalDofs=[1], dofValues=[P]) }

# Define the distributed loads (user-defined)
# for this problem, it is defined in frame's local coordinates, uniformly
# distributed.
def func_dload_loc(xloc):
    # define local x component of the applied traction function
    tx = 0 # no axial loading
    # define local y component of the applied traction function
    ty = -1 # N/mm # uniform vertical loading
    return [tx, ty]

# the following function defined on global x-y coords (for both inputs and outputs)
# produces the same result as the local one above. The key aspect is that this 
# function must be defined on the boundary itself, not along any projection of it;
# if it is the latter, then it must be transformed to a distribution along the 
# boundary while taking into account of the change of domain, i.e., the jacobian.
def func_dload_glb(xg):
    # set geometrical parameters
    L = 2000.0
    H = 1500.0
    # set ref point as the top_right point coords from the input file
    refP = [839.913107, 725.946989]
    # set tolerance for zero = 0.1 mm; it cannot be too small because it may 
    # exclude the gauss point, it just needs to be smaller than the distance 
    # between two gauss points and the distance between two elements
    tol = 1.e-1
    costheta = L/np.sqrt(L**2+H**2)
    sintheta = H/np.sqrt(L**2+H**2)
    q0 = 1
    tx = 0
    ty = 0
    # apply traction on the line
    if abs((xg[1]-refP[1])*L - (xg[0]-refP[0])*H) < tol:
        # define x component of the applied traction function
        tx = q0*sintheta
        # define y component of the applied traction function
        ty = -q0*costheta    
    return [tx, ty]

# 'local' indicates that tx and ty are along local frame axis and dload expression
# takes in point coordinates in local coordinate system
# order=0 indicates that the loading is constant
#dload_functions = [dload_function(expression=func_dload_loc, coord_system='local', order=0)]
dload_functions = [dload_function(expression=func_dload_glb, coord_system='global', order=0)]

#------------------------------------------------------------------------------
# define the interpreter connecting elset names to dload functions above 
#------------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! ADD/MODIFY THE SET NAMES BELOW TO BE CONSISTENT WITH INPUT FILE !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# all elements in the elset named 'q_elems' will be subjected to the distributed
# loading defined in dload_functions[0]
dict_elset_dload = {'FrameElems': dload_functions[0]}


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
# simple plotting within Python console for visualization of mesh
#------------------------------------------------------------------------------
r=100 # scaling factor for deformed plot
#plt.figure(figsize=[6.4,9.6])
plt.figure()
for elist in elem_lists:
    activeDofs = param.dict_eltype_activeDofs.get(elist.eltype)
    ndof_node  = len(activeDofs)
    for elem in elist.elems:
        elnodes = nodes[elem.cnc_node]
        eldofs  = a[elem.cnc_dof]
        # undeformed
        x_init = elnodes[:,0]
        y_init = elnodes[:,1]
        p1, =plt.plot(x_init,y_init,'.b-')
        # deformed
        x_curr = elnodes[:,0] + r*eldofs[0::ndof_node].reshape(len(elnodes))
        y_curr = elnodes[:,1] + r*eldofs[1::ndof_node].reshape(len(elnodes))
        p2, =plt.plot(x_curr, y_curr, '.r-')
p1.set_label('undeformed')
p2.set_label('deformed(scale=%d)' %r)
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.legend(loc='best')
plt.show()
#plt.savefig("deformedmeshscale100.png", dpi=250)
#plt.savefig("deformedmeshscale100.pdf")

# Derive strain & stress for elements
allstrain=[]
allstress=[]
for elist in elem_lists:
    if elist.eltype in param.tuple_truss2d2_eltypes:
        for elem in elist.elems:
            elnodes = nodes[elem.cnc_node]
            eldofs  = a[elem.cnc_dof].reshape(len(elnodes),2)
            eps, sig = truss2d2elem.strain_stress(elnodes[0], elnodes[1], \
                                eldofs[0], eldofs[1], Materials[elem.matID].E)
            allstrain.append([[eps]])
            allstress.append([[sig]])
            
    elif elist.eltype in param.tuple_frame2d2_eltypes:      
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
