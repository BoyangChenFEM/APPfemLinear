# -*- coding: utf-8 -*-
"""
Class for an integration point in a continuum element

@author: Boyang CHEN, TU Delft Aerospace, 2019
"""

class igpoint:
    
    def __init__(self, xi, w=None, x=None, u=None, strain=None, stress=None):
        self.xi = xi
        self.w  = w
        self.x  = x
        self.u  = u
        self.strain = strain
        self.stress = stress
        