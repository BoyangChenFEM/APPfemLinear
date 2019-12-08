# -*- coding: utf-8 -*-
"""
Created on Fri Nov 01 2019
Functions for a 2D 3-node linear triangular element
@author: Boyang CHEN TU Delft
"""
import numpy as np
import numpy.linalg as la

#class igpoint:
#    
#    def __init__(self, Nmatrix, strain, stress):
        

class tri2d3elem:
    
    def __init__(self, cnc_node, cnc_dof):
        self.cnc_node  = cnc_node
        self.cnc_dof   = cnc_dof
        
    @staticmethod
    def stiff_matrix(nodes, Dmat):
        """
        A function to calculate the stiffness matrix of the element 
        required inputs:
        nodes: coordinates of the nodes
        Dmat : material stiffness matrix
        """
        # calculate the area of the triangular element
        A = Area_triangle(nodes)
        # calculate the B matrix of the element
        B = B_matrix(nodes)
        # calculate the stiffness matrix K
        K = B.T@Dmat@B*A
        return K
    
    @staticmethod
    def fext_dload(nodes, dload_function):
        """
        TO BE DEFINED (PLACE HOLDER ONLY)
        A function to calculate the external force vector of the element due to 
        arbitrary distributed load along the edges of the element
        required inputs:
        nodes : global x-y coordinates of the nodes of the element
        dload_function : the distributed load vector as a function of (x,y)
        """

    
    @staticmethod
    def strain_stress(nodes, Dmat, a):
        """
        A function to calculate the strain and stress
        required inputs:
        nodes : coordinates of the nodes
        Dmat  : material stiffness matrix
        a     : dof vector of the element
        """
        # calculate the B matrix of the element
        B = B_matrix(nodes)
        # strain = B*a
        epsilon = B@a
        sigma   = Dmat@epsilon
        return [epsilon, sigma]



def Area_triangle(nodes):
    """
    A function to calculate the area for a triangular element
    input: nodal coordinates of the 3 nodes
    """
    x1 = nodes[0,0]
    y1 = nodes[0,1]
    x2 = nodes[1,0]
    y2 = nodes[1,1]
    x3 = nodes[2,0]
    y3 = nodes[2,1]
    A = 0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
    #A = 0.5*((x2*y3-x3*y2)+(y2-y3)*x1+(x3-x2)*y1)
    return A



def B_matrix(nodes):
    """ 
    A function to calculate the B matrix for the element
    inputs: nodal coordinates of the 3 nodes
    """
    A = Area_triangle(nodes)
    x1 = nodes[0,0]
    y1 = nodes[0,1]
    x2 = nodes[1,0]
    y2 = nodes[1,1]
    x3 = nodes[2,0]
    y3 = nodes[2,1]
    #a1 = x2*y3 - x3*y2
    #a2 = x3*y1 - x1*y3
    #a3 = x1*y2 - x2*y1
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2
    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1
    dN1dx = b1/(2*A) 
    dN1dy = c1/(2*A)
    dN2dx = b2/(2*A)
    dN2dy = c2/(2*A)
    dN3dx = b3/(2*A)
    dN3dy = c3/(2*A)
    B = np.array([[dN1dx,     0, dN2dx,     0, dN3dx,     0],\
                  [    0, dN1dy,     0, dN2dy,     0, dN3dy],\
                  [dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx]])
    return B

def fext_dload_uniform(node1, node2, t):
    """
    A function to calculate the external force vector of the element due to 
    uniformly distributed load along the edges of the element
    required inputs:
    node1&2: global x-y coordinates of the 2 nodes of the edge under loading
    t      : the uniformly distributed load vector
    """
    # calculate the length of the edge under loading
    L = la.norm(node2 - node1)
    # calculate the force vector due to uniform distributed loading
    f = 0.5*L*np.array([[t[0]],\
                        [t[1]],\
                        [t[0]],\
                        [t[1]]])
    return f