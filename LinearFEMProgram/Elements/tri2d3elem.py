# -*- coding: utf-8 -*-
"""
Class and functions for a 2D 3-node linear triangular element

@author: Boyang CHEN, TU Delft Aerospace, 2019
"""
import numpy as np
import numpy.linalg as la
from .integration_point import igpoint
        

class tri2d3elem:
    """class for tri2d3 element, with components:
        - cnc_node: nodal connectivity of this element
        - cnc_dof: dof connectivity of this element
        - matID: material ID of this element, default=0
        - igpoints: integration points; only 1 igpoint for this element type
        at centroid, for postprocessing only """ 
    def __init__(self, cnc_node, cnc_dof, matID=0):
        self.cnc_node  = cnc_node
        self.cnc_dof   = cnc_dof
        self.matID     = matID
        self.igpoints  = [igpoint(xi=[1/3, 1/3])]
        
    @staticmethod
    def stiff_matrix(nodes, Material):
        """
        A function to calculate the stiffness matrix of the element 
        required inputs:
        nodes: coordinates of the nodes
        Material: provides the material stiffness matrix
        """
        # calculate the area of the triangular element
        A = Area_triangle(nodes)
        # calculate the B matrix of the element
        B = B_matrix(nodes)
        # calculate the stiffness matrix K
        K = B.T@Material.Dmatrix@B*A
        return K
    
    @staticmethod
    def fext_dload(nodes, dload):
        """
        TO BE DEFINED (PLACE HOLDER ONLY)
        A function to calculate the external force vector of the element due to 
        arbitrary distributed load along the edges of the element. 
        Needed only if the load function is complicated that the equivalent
        nodal forces are not obvious to be derived and numerical integration 
        is needed.
        required inputs:
        nodes : global x-y coordinates of the nodes of the element
        dload.function : the distributed load vector as a function of coordinates
        dload.coord_system : the coordinate system of the function & its inputs
        dload.order: order of the function, to decide on the num. of igpoints
        """
        # initialize fext
        fext = np.zeros([6,1])
        
        # a list defining the local nodal connectivity of edges
        edges = [[0,1],[1,2],[2,0]]
        edges_dof = [[0,1,2,3],[2,3,4,5],[4,5,0,1]]
        
        if dload.order == 0:
            allxi = [0]
            allwt = [2]
        else:
            # 2-point Gauss scheme integrates up to 3rd order polynomial
            # shape function along edge is order 1, so dload function can be up to 
            # 2nd-order        
            allxi = [-1.0/np.sqrt(3), 1.0/np.sqrt(3)]
            allwt = [1.0, 1.0]
        
        if dload.order > 2:
            print('WARNING: the Gauss quadrature scheme in tri2d3elem may not\
                  be sufficient for this dload!')   
        
        for edge, edge_dof in zip(edges, edges_dof):
            for xi, wt in zip(allxi, allwt):
                # get the shape function values at xi in natural space 
                # same for all three edges
                # NOTE: shape functions of the tri2d3elem localized along the 
                # edge are nothing but linear truss shape functions.
                N1 = (1-xi)/2
                N2 = (1+xi)/2
                # use shape function to obtain gauss point coords in physical 
                # space along this edge: interpolate between end nodes of edge
                xgauss = N1*nodes[edge[0]] + N2*nodes[edge[1]]
                # evaluate the dload vector at the gauss point in physical space
                if dload.coord_system == 'global':
                    t = dload.expression(xgauss)
                else:
                    print('WARNING: dload in local coords is not yet supported \
                          in tri2d3elem!')
                # calculate the length of the edge under loading
                L = la.norm(nodes[edge[1]] - nodes[edge[0]])
                # calculate the force vector contribution of this gauss point
                # = jac * N^T * t
                ft = 0.5*L*np.array([[N1*t[0]],\
                                     [N1*t[1]],\
                                     [N2*t[0]],\
                                     [N2*t[1]]])
                # add its contribution to fext
                fext[edge_dof] += wt*ft
        
        return fext

    def update_igpoints(self, nodes, Material, a):
        """
        A function to update igpoints of self
        required inputs:
        nodes    : coordinates of the nodes
        Material : provides material stiffness matrix
        a        : dof vector of the element
        """
        # calculate & update position and displacement vectors
        self.igpoints[0].x = [1/3*np.sum(nodes[:,0]), 1/3*np.sum(nodes[:,1])]
        self.igpoints[0].u = [1/3*np.sum(a[0::2]),    1/3*np.sum(a[1::2])]
        # calculate strain and stress
        [epsilon, sigma] = strain_stress(nodes, Material.Dmatrix, a)
        # update to igpoint components
        self.igpoints[0].strain = epsilon
        self.igpoints[0].stress = sigma

    




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