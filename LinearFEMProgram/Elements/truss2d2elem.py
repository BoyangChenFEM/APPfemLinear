# -*- coding: utf-8 -*-
"""
Class and functions for a 2D 2-node truss element

@author: Boyang CHEN, TU Delft Aerospace, 2019
"""
import numpy as np
import numpy.linalg as la
from ..Materials import linear_elastic
        

class truss2d2elem:
    """class for truss2d2 element, with components:
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
        if not isinstance(Material, linear_elastic.truss):
            raise TypeError('Material type in truss2d2elem must be \
                            linear_elastic.truss')

    @staticmethod
    def stiff_matrix(nodes, Material):
        """ A function to calculate the stiffness matrix 
         required inputs:
         nodes : global x-y coordinates of the two end nodes
         Material: truss material which provides:
            - E : Young's modulus
            - A : cross-sectional area"""
        # extract materials constant
        E = Material.E
        A_cross = Material.A
        # calculate the length of the bar element
        L_e = la.norm(nodes[1] - nodes[0])
        # calculate the T matrix
        T = T_matrix(nodes[0], nodes[1])
        # calculate the local stiffness matrix K_e before rotation
        K_e = E*A_cross/L_e* \
        np.array([[ 1,  0, -1,  0],\
                  [ 0,  0,  0,  0],\
                  [-1,  0,  1,  0],\
                  [ 0,  0,  0,  0]])
        # rotate the local stiffness matrix to x-y coordinates
        K_e = T@K_e@T.T
        return K_e
    
    @staticmethod
    def fext_dload(nodes, dload):
        """
        nodes : global x-y coordinates of the nodes of the element
        dload : the distributed load data as object dload_data
        """
        fext = np.zeros([4,1])

        # calculate the length of the truss element
        L = la.norm(nodes[1] - nodes[0])
        # jacobian from physical to natural space
        jac = L/2
        
        if dload.order == 0:
            allxi = [0]
            allwt = [2]
        else:
            # 2-point Gauss scheme integrates up to 3rd order polynomial
            # shape function along truss is order 1, so dload function can be up to 
            # 2nd-order        
            allxi = [-1.0/np.sqrt(3), 1.0/np.sqrt(3)]
            allwt = [1.0, 1.0]
        
        if dload.order > 2:
            print('WARNING: the Gauss quadrature scheme in truss2d2elem may not\
                  be sufficient for this dload!') 
            
        for xi, wt in zip(allxi, allwt):
            # get the shape function values at xi in natural space
            N1 = (1-xi)/2
            N2 = (1+xi)/2
            # use shape function to obtain gauss point coords in physical 
            # space between end nodes
            xg = N1*nodes[0] + N2*nodes[1]
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
                ft = T@NTq(N1, N2, q)
            elif dload.coord_system == 'global':
                q = dload.expression(xg)
                # calculate the force vector contribution of this gauss point
                ft = NTq(N1, N2, q)
            else:
                print('WARNING: dload coord system is not supported \
                      in truss2d2elem!')
            
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
    # fill in the T matrix
    T = np.array([[costheta,  -sintheta,  0,  0],\
                  [sintheta,   costheta,  0,  0],\
                  [0,  0,  costheta,  -sintheta],\
                  [0,  0,  sintheta,   costheta]])
    return T


def NTq(N1, N2, q):
    # calculate the force vector contribution N^T * q
    ft = np.array([[N1*q[0]],\
                   [N1*q[1]],\
                   [N2*q[0]],\
                   [N2*q[1]]])
    return ft


def strain_stress(node1, node2, u1, u2, E):
    """A function to calculate the strain and stress 
    required inputs:
    node 1 & 2 : global x-y coordinates of the two end nodes
    u1 & 2     : displacement vectors of the two end nodes in x-y coordinates
    E          : Young's modulus"""
    # calculate the initial length of the bar element
    L_0 = la.norm(node2 - node1)
    # calculate the current length of the bar element
    L_t = la.norm(node2 + u2 - node1 - u1)
    # calculate the strain
    epsilon = (L_t-L_0)/L_0
    # calculate the stress
    sigma = E*epsilon
    return [epsilon,sigma]