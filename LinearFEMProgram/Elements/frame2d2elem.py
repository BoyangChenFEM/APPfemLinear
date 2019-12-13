# -*- coding: utf-8 -*-
"""
Created on Fri Nov 01 2019
Functions for a 2D 3-node linear triangular element
@author: Boyang CHEN TU Delft
"""
import numpy as np
import numpy.linalg as la
        

class frame2d2elem:
    
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
        Material : provides material stiffness matrix/constants
        """
        E = Material.E
        A_cross = Material.A
        I = Material.I
        # calculate the length of the bar element
        L = la.norm(nodes[1] - nodes[0])
        # calculate the transformation matrix
        T = T_matrix(nodes[0], nodes[1])
        # calculate the local stiffness matrix K_e before rotation
        # K_t: truss part of the stiffness matrix
        # K_b: beam part of the stiffness matrix
        K_t = E*A_cross/L* \
        np.array([[ 1,  0,  0, -1,  0,  0],\
                  [ 0,  0,  0,  0,  0,  0],\
                  [ 0,  0,  0,  0,  0,  0],\
                  [-1,  0,  0,  1,  0,  0],\
                  [ 0,  0,  0,  0,  0,  0],\
                  [ 0,  0,  0,  0,  0,  0]])
        K_b = E*I* \
        np.array([[ 0,        0,       0,  0,        0,       0],\
                  [ 0,  12/L**3,  6/L**2,  0, -12/L**3,  6/L**2],\
                  [ 0,   6/L**2,     4/L,  0,  -6/L**2,     2/L],\
                  [ 0,        0,       0,  0,        0,       0],\
                  [ 0, -12/L**3, -6/L**2,  0,  12/L**3, -6/L**2],\
                  [ 0,   6/L**2,     2/L,  0,  -6/L**2,     4/L]])
        # rotate the total local stiffness matrix to x-y coordinates
        K = T@(K_t+K_b)@T.T
        return K
    
    @staticmethod
    def fext_dload(nodes, dload):
        """
        nodes : global x-y coordinates of the nodes of the element
        dload : the distributed load vector as a function of (x,y)
        """
        fext = np.zeros([6,1])
        if dload.coord_system == 'local' and dload.order == 0:
            fext = fext_dload_uniform(nodes, dload.expression)
        else:
            print('WARNING: global or non-uniformly distributed load on frames\
                  is not yet implemented.')
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



def strain_stress(nodes, a_glb, E, h, point):
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
    # Bending B matrix: Bb = d^2 N_bending / dx^2
    Bb = np.array([0, 12*x/L**3-6/L**2, 6*x/L**2-4/L, 0, 6/L**2-12*x/L**3, 6*x/L**2-2/L])
    # Axial strain = Ba*a, Bending strain = - y*Bb*a
    epsilon = Ba@a - y*Bb@a
    sigma   = E*epsilon
    return [epsilon,sigma]