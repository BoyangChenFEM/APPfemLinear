# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 22:05:38 2019

@author: Boyang
"""

class igpoint:
    
    def __init__(self, xi, w, x=None, u=None, strain=None, stress=None):
        self.xi = xi
        self.w  = w
        self.x  = x
        self.u  = u
        self.strain = strain
        self.stress = stress
        