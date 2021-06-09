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
    def verify_material(Material):
        """ check if the material type is correct """
        if not isinstance(Material, linear_elastic.frame2D):
            raise TypeError('Material type in frame2d2elem must be \
                            linear_elastic.frame2D')
        
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
        
        # define the stiffness matrix   
        # first, calculate the truss part
        K_t = K_truss(Material.E, Material.A, L)
        # next, calculate the beam part
        K_b = K_beam(Material.E, Material.I, L)
        # rotate the total local stiffness matrix to x-y coordinates
        K = T@(K_t+K_b)@T.T
            
        return K
    
    @staticmethod
    def fext_dload(nodes, dload):
        """
        nodes : global x-y coordinates of the nodes of the element
        dload : the distributed load data as object dload_data
        """
        fext = np.zeros([6,1])
        
        # calculate the length of the frame element
        L = la.norm(nodes[1] - nodes[0])
        # jacobian from physical to natural space
        jac = L/2
        
        if dload.order == 0:
            # 2-point Gauss scheme integrates up to 3rd order polynomial
            # shape function along frame is order 3, so dload function can be up to 
            # 0th-order        
            allxi = [-1/np.sqrt(3), 1/np.sqrt(3)]
            allwt = [1, 1]
        else:
            # 3-point Gauss scheme integrates up to 5th order polynomial
            # shape function along frame is order 3, so dload function can be up to 
            # 2nd-order        
            allxi = [-np.sqrt(0.6), 0, np.sqrt(0.6)]
            allwt = [5.0/9.0, 8.0/9.0, 5.0/9.0]
        
        if dload.order > 2:
            print('WARNING: the Gauss quadrature scheme in frame2d2elem may not\
                  be sufficient for this dload!')
        
        for xi, wt in zip(allxi, allwt):
            # get the truss shape function values at xi in natural space 
            N1t = (1-xi)/2
            N2t = (1+xi)/2
            # use truss shape function to obtain gauss point coords in physical 
            # space along this frame: interpolate between end nodes of frame
            xg = N1t*nodes[0] + N2t*nodes[1]
            # calculate local coord xhat of gauss point along frame within [0, L]
            xhat = N2t*L
            # calculate the shape function values of beam part
            N1b = 1 - 3*(xhat/L)**2 + 2*(xhat/L)**3
            N2b = xhat - 2*xhat**2/L + xhat**3/L**2
            N3b = 3*(xhat/L)**2 - 2*(xhat/L)**3
            N4b = xhat**3/L**2 - xhat**2/L
#            # calculate the N matrix of the frame element at this gauss point
#            Nmat = np.array([ [N1t,  0,   0, N2t,   0,   0],\
#                              [0,  N1b, N2b,   0, N3b, N4b] ])
            # evaluate the dload vector at the gauss point in physical space
            if dload.coord_system == 'local':
                # calculate the transformation matrix Q
                Q = Q_matrix(nodes[0], nodes[1])
                # transform gauss point vector from global coords to local
                xl = Q.T@xg
                # pass local gauss point coords to load function expression
                # to obtain the traction vector in local coordinates
                q = dload.expression(xl)
                # calculate the force vector contribution in global coords
                T = T_matrix(nodes[0], nodes[1])
                ft = T@NTq(N1t, N2t, N1b, N2b, N3b, N4b, q)
            elif dload.coord_system == 'global':
                Q = Q_matrix(nodes[0], nodes[1])
                q = Q.T@dload.expression(xg)
                T = T_matrix(nodes[0], nodes[1])
                ft = T@NTq(N1t, N2t, N1b, N2b, N3b, N4b, q)
            else:
                print('WARNING: dload coord system is not supported \
                      in frame2d2elem!')

            # add its contribution to fext
            fext += jac*wt*ft
        
        return fext



def Q_matrix(node1, node2):
    """ A function to calculate the transformation matrix Q for the element.
    Q transforms a position vector from local to global coords: 
    r_glb = Q * r_lcl"""
    # calculate the unit axial vector of the bar element in global coords
    axial_vect = (node2 - node1)/la.norm(node2 - node1)
    # get cos theta and sin theta
    costheta = axial_vect[0]
    sintheta = axial_vect[1]
    # fill in the transformation matrix
    Q = np.array([[costheta, -sintheta],\
                  [sintheta,  costheta]])
    return Q
    

def T_matrix(node1, node2):
    """ A function to calculate the transformation matrix T for the element.
    T transforms a dof/force vector from local to global coords: 
    v_glb = T * v_lcl"""
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


def NTq(N1t, N2t, N1b, N2b, N3b, N4b, q):
    # calculate the force vector contribution N^T * q
    ft = np.array([[N1t*q[0]],\
                   [N1b*q[1]],\
                   [N2b*q[1]],\
                   [N2t*q[0]],\
                   [N3b*q[1]],\
                   [N4b*q[1]]])
    return ft


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
     point      : normalized local coordinate of the probing point, i.e., 
                 point[0] is between [0,1] and point[1] is between [-1/2,1/2]
                 for beams with neutral axis in the center.
     """
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
    sigma   = Material.E*epsilon
    return [epsilon,sigma]