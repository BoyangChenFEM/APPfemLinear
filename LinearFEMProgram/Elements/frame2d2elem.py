# -*- coding: utf-8 -*-
"""
Class and functions for a 2D 2-node frame element

@author: Boyang CHEN, TU Delft Aerospace, 2019
"""
import numpy as np
import numpy.linalg as la
from ..Materials import linear_elastic
        

class frame2d2elem:
    """class for frame2d2 element, with components:
        - cnc_node: nodal connectivity of this element
        - cnc_dof: dof connectivity of this element
        - matID: material ID of this element, default=0"""    
    def __init__(self, cnc_node, cnc_dof, matID=0):
        self.cnc_node  = cnc_node
        self.cnc_dof   = cnc_dof
        self.matID     = matID
        
    @staticmethod
    def stiff_matrix(nodes, Material):
        """
        A function to calculate the stiffness matrix of the element 
        required inputs:
        nodes: coordinates of the nodes
        Material : provides material stiffness matrix/constants; it can be both
        truss and frame in type
        """
        # calculate the length of the bar element
        L = la.norm(nodes[1] - nodes[0])
        # calculate the transformation matrix
        T = T_matrix(nodes[0], nodes[1])
        
        # define the stiffness matrix based on material type
        # here, both truss and frame2D materials are supported
        
        # first, calculate the truss part in any case
        K_t = K_truss(Material.E, Material.A, L)
        
        # form K based on material type
        if isinstance(Material, linear_elastic.truss):
            # add dummy stiffness 1 to the diagonal terms on theta1 & theta2
            # to avoid singular matrix
            K_t[2,2] = 1
            K_t[5,5] = 1
            # rotate the total local stiffness matrix to x-y coordinates
            K = T@K_t@T.T
        elif isinstance(Material, linear_elastic.frame2D):
            K_b = K_beam(Material.E, Material.I, L)
            # rotate the total local stiffness matrix to x-y coordinates
            K = T@(K_t+K_b)@T.T
        else:
            print('WARNING: unsupported material in frame2delem')
            
        return K
    
    @staticmethod
    def fext_dload(nodes, dload):
        """
        nodes : global x-y coordinates of the nodes of the element
        dload : the distributed load data as object dload_data
        """
        fext = np.zeros([6,1])
        
#        # 3-point Gauss scheme integrates up to 5th order polynomial
#        # shape function along frame is order 3, so dload function can be up to 
#        # 2nd-order        
#        allxi = [-np.sqrt(0.6), 0, np.sqrt(0.6)]
#        allwt = [5.0/9.0, 8.0/9.0, 5.0/9.0]
#        if dload.order > 2:
#            print('WARNING: the Gauss quadrature scheme in frame2d2elem may not\
#                  be sufficient for this dload!')
#        
#        for xi, wt in zip(allxi, allwt):
#            # get the truss shape function values at xi in natural space 
#            N1 = (1-xi)/2
#            N2 = (1+xi)/2
#            # use truss shape function to obtain gauss point coords in physical 
#            # space along this frame: interpolate between end nodes of frame
#            xgauss = N1*nodes[0] + N2*nodes[1]
#            # evaluate the dload vector at the gauss point in physical space
#            if dload.coord_system == 'global':
#                q = dload.expression(xgauss)
#                # TO BE CONTINUED
#                
#            elif dload.coord_system == 'local':
#                # calculate the transformation matrix
#                T = T_matrix(nodes[0], nodes[1])
#                # transform gauss point from global coords to local
#                xgauss = T.T@xgauss
#                # pass local gauss point coords to load function expression
#                q = dload.expression(xgauss)
#                # TO BE CONTINUED
#                
#            else:
#                print('WARNING: dload coord system is not supported \
#                      in frame2d2elem!')
##            # calculate the length of the edge under loading
##            L = la.norm(nodes[edge[1]] - nodes[edge[0]])
##            # calculate the force vector contribution of this gauss point
##            # = jac * N^T * t
##            ft = 0.5*L*np.array([[N1*t[0]],\
##                                 [N1*t[1]],\
##                                 [N2*t[0]],\
##                                 [N2*t[1]]])
#            # add its contribution to fext
#            fext += wt*ft
        
#        if dload.coord_system == 'local' and dload.order == 0:
#            fext = fext_dload_uniform(nodes, dload.expression())
#        else:
#            print('WARNING: global or non-uniformly distributed load on frames\
#                  is not yet implemented.')
        return fext


    


def T_matrix(node1, node2):
    """ A function to calculate the transformation matrix T for the element
    T transforms a vector from local to global coords: v_glb = T * v_lcl"""
    # calculate the unit axial vector of the bar element in global coords
    axial_vect = (node2 - node1)/la.norm(node2 - node1)
    # get cos theta and sin theta
    costheta = axial_vect[0]
    sintheta = axial_vect[1]
    # fill in the transformation matrix
    T = np.array([[costheta, -sintheta,  0,         0,          0,  0],\
                  [sintheta,  costheta,  0,         0,          0,  0],\
                  [       0,         0,  1,         0,          0,  0],\
                  [       0,         0,  0,  costheta,  -sintheta,  0],\
                  [       0,         0,  0,  sintheta,   costheta,  0],\
                  [       0,         0,  0,         0,          0,  1]])
    return T


def K_truss(E, A_cross, L):
    K_t = E*A_cross/L* \
    np.array([[ 1,  0,  0, -1,  0,  0],\
              [ 0,  0,  0,  0,  0,  0],\
              [ 0,  0,  0,  0,  0,  0],\
              [-1,  0,  0,  1,  0,  0],\
              [ 0,  0,  0,  0,  0,  0],\
              [ 0,  0,  0,  0,  0,  0]])
    return K_t


def K_beam(E, I, L):
    K_b = E*I* \
    np.array([[ 0,        0,       0,  0,        0,       0],\
              [ 0,  12/L**3,  6/L**2,  0, -12/L**3,  6/L**2],\
              [ 0,   6/L**2,     4/L,  0,  -6/L**2,     2/L],\
              [ 0,        0,       0,  0,        0,       0],\
              [ 0, -12/L**3, -6/L**2,  0,  12/L**3, -6/L**2],\
              [ 0,   6/L**2,     2/L,  0,  -6/L**2,     4/L]])
    return K_b


def fext_dload_uniform(nodes, q):
    """ A function to calculate the external force vector of the element due to 
     uniformly distributed load defined in local coordinates
     required inputs:
     node 1 & 2 : global x-y coordinates of the two end nodes
     q          : vector of the uniformly distributed load"""
    # calculate the length of the bar element
    L = la.norm(nodes[1] - nodes[0])
    # calculate the transformation matrix
    T = T_matrix(nodes[0], nodes[1])
    # calculate the local force vector due to uniform pressure
    f_lcl = np.array([
              [q[0]*L/2],\
              [q[1]*L/2],\
              [q[1]*L**2/12],\
              [q[0]*L/2],\
              [q[1]*L/2],\
              [-q[1]*L**2/12]
              ])
    # rotate the total local stiffness matrix to x-y coordinates
    f = T@f_lcl
    return f



def strain_stress(nodes, a_glb, Material, h, point):
    """ A function to calculate the strain and stress
     required inputs:
     node 1 & 2 : global x-y coordinates of the two end nodes
     a_glb      : global dof vector of the element in x-y coordinates
     E          : Young's modulus
     point      : local coordinate of the probing point"""
    # calculate the length of the bar element
    L = la.norm(nodes[1] - nodes[0])
    # calculate the transformation matrix
    T = T_matrix(nodes[0], nodes[1])
    # local dof vector, v_lcl = transpose(T) * v_glb
    a = T.T@a_glb
	# local coordinates of the point of stress/strain evaluation
    x = point[0]*L
    y = point[1]*h
    # axial B matrix, using the gradient of the axial shape function 
    Ba = np.array([-1/L, 0, 0, 1/L, 0, 0])
    
    if isinstance(Material, linear_elastic.truss):
        epsilon = Ba@a
    elif isinstance(Material, linear_elastic.frame2D):
        # Bending B matrix: Bb = d^2 N_bending / dx^2
        Bb = np.array([0, 12*x/L**3-6/L**2, 6*x/L**2-4/L, 0, 6/L**2-12*x/L**3, 6*x/L**2-2/L])
        # Axial strain = Ba*a, Bending strain = - y*Bb*a
        epsilon = Ba@a - y*Bb@a
    else:
        print('WARNING: unsupported material in frame2delem')
    
    sigma   = Material.E*epsilon
    return [epsilon,sigma]