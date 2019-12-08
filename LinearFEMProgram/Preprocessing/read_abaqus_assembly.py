# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 17:26:43 2019

Abaqus input file reader for assembly

@author: Boyang CHEN TU Delft
"""
from read_abaqus_parts import read_set_name, read_set_list

class assembly:
    
    def __init__(self, name, instances, nsets=[], elsets=[]):
        self.name = name
        self.instances = instances
        self.nsets  = nsets
        self.elsets = elsets

class instance:

    def __init__(self, name, part):
        self.name = name
        self.part = part        
    
class assembly_set:

    def __init__(self, name, instance, setlist):
        self.name  = name
        self.instance = instance
        self.nodes = setlist

        
# =================================================
# Define functions to read objects in parts from
# Abaqus input file
# =================================================

def form_assembly(lines):
    instances = []
    nsets     = []
    elsets    = []
    # find the line no. of *Assembly and *End Assembly and store them
    jstart = next( j for j,line in enumerate(lines) if '*Assembly' in line )
    jend   = next( j for j,line in enumerate(lines) if '*End Assembly' in line )
    # read assembly name
    name = lines[jstart].rstrip().split(',')[1].split('=')[1]
    # read instances
    jinstances = [ j for j,line in enumerate(lines[jstart:jend]) if '*Instance' in line ]
    for jinstance in jinstances:
        instances.append(read_instance(lines[jinstance]))
    # read nsets
    jnsets = [ j for j,line in enumerate(lines[jstart:jend]) if '*Nset' in line ]
    for jnset in jnsets:
        nsets.append(read_assembly_set(lines[jnset:jend]))
    # read elsets
    jelsets = [ j for j,line in enumerate(lines[jstart:jend]) if '*Elset' in line ]
    for jelset in jelsets:
        elsets.append(read_assembly_set(lines[jelset:jend]))
    # form the assembly and return    
    return assembly(name, instances, nsets, elsets)


def read_instance(line): 
    name = line.rstrip().split(',')[1].split('=')[1]
    part = line.rstrip().split(',')[2].split('=')[1]
    return instance(name, part)


def read_assembly_set(lines):
    # get the set name
    line0    = lines[0].rstrip()
    name     = read_set_name(line0)  
    instance = read_set_instance(line0)
    setlist  = read_set_list(lines)
    return assembly_set(name, instance, setlist)


def read_set_instance(line):
    instance = line.split(',')[2].split('=')[1]
    return instance