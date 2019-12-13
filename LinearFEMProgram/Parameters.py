# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 10:39:08 2019
A module to store problem-independent parameters. 
They should be used by the Kernel only.
@author: Boyang
"""
# no. of decimal places for output purposes
num_decimal_for_output = 6

# A dictionary to relate Abaqus element type and its no. of dofs per node
dict_eltype_ndofnode = {'T2D2':2, 'B23':3, \
                        'CPS3':2, 'CPE3':2, \
        'CPS4':2, 'CPS4R':2, 'CPE4':2, 'CPE4R':2}

# A dictionary to relate Abaqus element type and its no. of dimensions
dict_eltype_ndim = {'T2D2':2, 'B23':2, \
                    'CPS3':2, 'CPE3':2, \
        'CPS4':2, 'CPS4R':2, 'CPE4':2, 'CPE4R':2}

# A dictionary to relate Abaqus element type and its no. of strains
dict_eltype_nst = {'T2D2':1, 'B23':1, \
                   'CPS3':3, 'CPE3':3, \
        'CPS4':3, 'CPS4R':3, 'CPE4':3, 'CPE4R':3}

# tuples of supported element types
tuple_truss2d2_eltypes = ('T2D2',)
tuple_frame2d2_eltypes = ('B23',)
tuple_tri2d3_eltypes = ('CPS3', 'CPE3')
#tuple_quad2d4_eltypes = ('CPS4', 'CPS4R', 'CPE4', 'CPE4R')
tuple_supported_eltypes = tuple_truss2d2_eltypes + \
                          tuple_frame2d2_eltypes + \
                          tuple_tri2d3_eltypes

