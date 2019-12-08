# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 23:40:59 2019

@author: Boyang
"""
import numpy as np

def Dmatrix_2D_isotropic(E, nu, is_planestress):
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
    return Dmat