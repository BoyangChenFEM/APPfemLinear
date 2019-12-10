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
dict_eltype_ndofnode = {'CPS3':2, 'CPE3':2, \
        'CPS4':2, 'CPS4R':2, 'CPE4':2, 'CPE4R':2}

# A dictionary to relate Abaqus element type and its no. of dimensions
dict_eltype_ndim = {'CPS3':2, 'CPE3':2, \
        'CPS4':2, 'CPS4R':2, 'CPE4':2, 'CPE4R':2}

# A dictionary to relate Abaqus element type and its no. of strains
dict_eltype_nst = {'CPS3':3, 'CPE3':3, \
        'CPS4':3, 'CPS4R':3, 'CPE4':3, 'CPE4R':3}

# tuples of supported element types
tuple_tri2d3_eltypes = ('CPS3', 'CPE3')
#tuple_quad2d4_eltypes = ('CPS4', 'CPS4R', 'CPE4', 'CPE4R')
tuple_supported_eltypes = tuple_tri2d3_eltypes #+ tuple_quad2d4_eltypes

