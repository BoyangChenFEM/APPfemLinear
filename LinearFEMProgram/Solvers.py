# -*- coding: utf-8 -*-
"""
Linear FEM Assembler and Solver

They are put in the same module because they need to be consistent in 
the data format of K (sparse or not)

@author: Boyang CHEN, TU Delft Aerospace, 2019
"""
import numpy as np
from scipy.sparse import lil_matrix # sparse matrix format as linked list
import scipy.linalg as la
from pypardiso import spsolve #Intel MKL Pardiso parrallel(shared memory) sparse matrix solver

def assembler(nodes, elems, NDOF_NODE, Materials, list_dload_data=[]):
    """Assemble the system equation K a = f, where K is assembly of element 
    stiffness matrix and f is assembly of the element external force vector 
    due to distributed loading (if present); note that the concentrated forces
    will be applied directly on f in the solver """
    # calculate the size of the system
    # total no. of dofs in the system
    ndofs = NDOF_NODE*len(nodes)
    
    # Initialize stiffness matrix K and force vector
    #K = np.zeros([ndofs,ndofs])
    K = lil_matrix((ndofs,ndofs)) #sparse matrix, lil format
    f = np.zeros([ndofs,1])
    
    # Calculate element stiffness matrix K_e and assemble to K
    for elem in elems:    
        # extract the nodal coordinates of the element using cnc matrix
        elem_nodes = nodes[elem.cnc_node]
        # calculate the element stiffness matrix
        K_e = elem.stiff_matrix(elem_nodes, Materials[elem.matID])
        # assemble K_e to K using cnc_dof matrix
        K[np.ix_(elem.cnc_dof[:],elem.cnc_dof[:])] += K_e[:,:]
    
    # get f_ext due to distributed loading
    if list_dload_data:
        for dload in list_dload_data:
            for ie in dload.elset:
                elem_nodes = nodes[elems[ie].cnc_node]
                f_e = elems[ie].fext_dload(elem_nodes, dload.function)
                f[elems[ie].cnc_dof] += f_e
    
    return [K, f]


def direct_solver(K, f, bcd_dofs, bcd_values, cload_dofs, cload_values):
    """Direct solver of fem system equation by decomposition of dof vector into
    free part and constraint part in the style of Tom Hughes"""
    # get ndofs: total no. of dofs
    ndofs = len(f)
    
    # initialize solution vector a
    a = np.zeros([ndofs,1])

    # get free dofs, i.e., dofs without imposed bcds
    free_dofs = np.setdiff1d(list(range(ndofs)),bcd_dofs)    

    # apply the equivalent force vector on free dofs due to imposed displacements 
    f[free_dofs,0] -= K[np.ix_(free_dofs,bcd_dofs)]@np.asarray(bcd_values)    

    # Apply the concentrated loads on the force vector
    f[cload_dofs,0] += cload_values
    
    # identify zero rows in K and put diagonal to 1 to avoid singularity
    #zeroRows = np.where(~K.any(axis=1))[0]
    #K[zeroRows, zeroRows] = 1 # set zero row&col diagonal term to 1    
    zeroRows = np.diff(K.tocsr().indptr)==0
    K[zeroRows, zeroRows] = 1 # set zero row&col diagonal term to 1
    
    # Solve for a:
    #a[free_dofs] = la.solve(K[np.ix_(free_dofs,free_dofs)], f[free_dofs])
    a[free_dofs,0] = spsolve(K[np.ix_(free_dofs,free_dofs)].tocsr(), f[free_dofs])
	
    # update a to include the imposed displacements
    a[bcd_dofs,0] = bcd_values
    
    # obtain reaction force at node with imposed displacement
    RF = np.zeros([ndofs,1])
    RF[bcd_dofs] = K[bcd_dofs,:]@a
    
    return [a, RF]



def direct_solver_edu(K, f, bcd_dofs, bcd_values, cload_dofs, cload_values):
    """Educational version of the direct solver, showing crossing of rows and 
    columns for direct connections with fem instructions in the style of 
    Zienkiewwicz & Zhu, Liu & Quek"""
    # get ndofs: total no. of dofs
    ndofs = len(f)
    
    # cross-out rows with imposed displacement; here set them all to zero
    # store this row for postprocessing later before crossing it out
    Krow_bcd = np.zeros([len(bcd_dofs),ndofs])
    Krow_bcd = K[bcd_dofs[:],:] # storage
    K[bcd_dofs[:],:] = 0 # cross-out row in K
    
    # if imposed displacement is not zero, then move U_D * column(bcd_node) to the
    # RHS and set the column wrt U_D to zero
    for x,i in zip(bcd_values,bcd_dofs):
        if (x != 0):
            f[:] -= x * K[:,i]
    
    # cross-out column of K wrt node with imposed displacement
    K[:,bcd_dofs[:]] = 0
    
    # Apply the concentrated loads on the force vector
    f[cload_dofs[:],0] += cload_values[:]
    
#    # set diagonal term of K wrt bcd node to a dummy stiffness to avoid singularity
#    for i in bcd_dofs:
#       K[i,i] = 1 # dummy non-zero stiffness 

    # identify zero rows in K and put diagonal to 1 to avoid singularity
    # this operation makes the above-commented section redundant
    #zeroRows = np.where(~K.any(axis=1))[0]
    #K[zeroRows, zeroRows] = 1 # set zero row&col diagonal term to 1    
    zeroRows = np.diff(K.tocsr().indptr)==0
    K[zeroRows, zeroRows] = 1 # set zero row&col diagonal term to 1
    
    # Solve for a:
    #a = la.solve(K,f)
    a = spsolve(K.tocsr(),f).reshape(ndofs,1)
    
    # update a to include the imposed displacements
    a[bcd_dofs[:],0] = bcd_values[:]
    
    # obtain reaction force at nodes with imposed displacements
    RF = np.zeros([ndofs,1])
    RF[bcd_dofs] = Krow_bcd@a
    
    return [a, RF]