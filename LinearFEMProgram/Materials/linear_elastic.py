# -*- coding: utf-8 -*-
"""
A module of classes and methods for linear elastic materials and sections

@author: Boyang CHEN, TU Delft Aerospace, 2019
"""
import numpy as np

class truss:
    def __init__(self, E, A):
        self.E = E
        self.A = A
        
class beam2D:
    def __init__(self, E, I):
        self.E = E
        self.I = I

class frame2D:
    def __init__(self, E, A, I):
        self.E = E
        self.A = A
        self.I = I

class isotropic2D:
    def __init__(self, E, nu, thickness, is_planestress):
        self.E = E
        self.nu = nu
        self.thickness = thickness
        self.is_planestress = is_planestress
        self.Dmatrix = Dmatrix_2D_isotropic(E, nu, thickness, is_planestress)

def Dmatrix_2D_isotropic(E, nu, thickness, is_planestress):
    """A function to get the 2D linear elastic isotropic material stiffness
    matrix"""
    if is_planestress:
        Dmat = E/(1-nu**2)*\
               np.array([[ 1, nu,         0],\
                         [nu,  1,         0],\
                         [ 0,  0,  (1-nu)/2]])
    else:
        Dmat = E/((1+nu)*(1-2*nu))*\
               np.array([[1-nu,    nu,           0],\
                         [  nu,  1-nu,           0],\
                         [   0,     0,  (1-2*nu)/2]])
    return thickness*Dmat