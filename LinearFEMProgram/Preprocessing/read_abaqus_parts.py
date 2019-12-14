# -*- coding: utf-8 -*-
"""
Abaqus input file reader for parts

@author: Boyang CHEN, TU Delft Aerospace, 2019
"""

# =================================================
# Define classes for objects in abaqus parts
# =================================================

class elem_group:
    
    def __init__(self, eltype, cnc_node):
        self.eltype    = eltype
        self.cnc_node  = cnc_node
   
class part_set:

    def __init__(self, name, setlist):
        self.name    = name
        self.setlist = setlist

class part:

    def __init__(self, name, nodes, elem_groups, nsets=[], elsets=[]):
        self.name        = name
        self.nodes       = nodes
        self.elem_groups = elem_groups
        self.nsets       = nsets
        self.elsets      = elsets


# =================================================
# Define functions to read objects in parts from
# Abaqus input file
# =================================================

def read_parts_from_inputfile(inputfile):
    # Read data from Abaqus input file and form abaqus parts
    with open(inputfile,'r') as abaqus_input:
        # Store all lines first
        lines = abaqus_input.readlines()
        # define parts; parts is an object defined in read_abaqus_parts module which 
        # stores the *Part information in the abaqus input file
        parts = form_parts(lines)
        # close abaqus input file
        abaqus_input.close()
    return parts


def form_parts(lines):
    # initialize parts
    parts = []
    # find the line no. of *Part and *End Part and store them
    jparts_start = [j for j,line in enumerate(lines) if '*part' in line.lower()]
    jparts_end   = [j for j,line in enumerate(lines) if '*end part' in line.lower()]
    # if no part or no end part, then form one part only
    if (not jparts_start or not jparts_end):
        parts.append(form_a_part(lines))
    # read all parts
    else:
        for jps,jpe in zip(jparts_start,jparts_end):
            # read one part and append to the parts list
            parts.append(form_a_part(lines[jps:jpe]))
    return parts


def form_a_part(lines):
    # read Part name
    if '*part' in lines[0].lower():
        pname = lines[0].rstrip().split(',')[1].split('=')[1]
    else:
        pname = 'Part-1'
    # find the line of *Node, *Element, *Nset, *Elset
    # only ONE Node section (default)
    jnodes = next( j for j in range(len(lines))  if '*node'      in lines[j].lower() )
    jelems =     [ j for j in range(len(lines))  if '*element'   in lines[j].lower() ]
    jnsets =     [ j for j in range(len(lines))  if '*nset'      in lines[j].lower() ]
    jlsets =     [ j for j in range(len(lines))  if '*elset'     in lines[j].lower() ]
    # read nodes
    nodes = read_nodes(lines[jnodes+1:])
    # read elements
    elem_groups = []
    for je in jelems:
        elem_groups.append(read_elems(lines[je:]))
    # read nsets
    nsets = []
    for jnset in jnsets:
        nsets.append(read_part_set(lines[jnset:]))
    # read elsets
    elsets = []
    for jlset in jlsets:
        elsets.append(read_part_set(lines[jlset:]))    
    # define and return part based on part components
    return part(pname, nodes, elem_groups, nsets, elsets)


def read_nodes(lines):
    nodes = []
    for line in lines:
        # break out of for loop if end of node section is reached
        if('*' in line):
            break
        # read the coords of this node into a list of float numbers
        coords = []
        for t in line.split(','):
            try:
                coords.append(float(t))
            except ValueError:
                pass
        # store the coords in nodes list of this part (coords[0] is the node no.)
        nodes.append(coords[1:])
    return nodes


def read_elems(lines):
    # read element type
    eltype = lines[0].rstrip().split(',')[1].split('=')[1]
    # read element connectivities
    cnc_node = []
    for line in lines[1:]:
        # break out of for loop if end of elem section is reached
        if('*' in line):
            break
        # read the index and nodes of this elem into a list of int numbers
        el = []
        for t in line.split(','):
            try:
                el.append(int(t)-1) # -1 to start nunbering from 0
            except ValueError:
                pass
        # el[0]  : elem index
        # el[1:] : elem nodes
        cnc_node.append(el[1:])
    return elem_group(eltype, cnc_node)


def read_part_set(lines):
    # get the set name
    line0 = lines[0].rstrip()
    set_name = read_set_name(line0)
    set_list = read_set_list(lines)
    return part_set(set_name, set_list)


def read_set_name(line):
    set_name = line.split(',')[1].split('=')[1]
    return set_name


def read_set_list(lines):
    # create empty list
    set_list = []
    # read set
    # if generate is used, then calculate all nodes;
    # otherwise, read all nodes directly
    if ('generate' in lines[0].lower()):
        nl = [] # node list
        line = lines[1]
        for t in line.split(','):
            try:
                nl.append(int(t)) 
            except ValueError:
                pass
        nds = nl[0]-1 # start node # -1 to start nunbering from 0 
        ndf = nl[1]-1 # final node # -1 to start nunbering from 0 
        try:
            itv = nl[2] # interval
        except IndexError:
            itv = 1 # default interval is 1
        for n in range(nds,ndf+1,itv):
            set_list.append(n) 
    else:
        nl = [] # node list
        for line in lines[1:]:
            # break out of loop if end of section encountered
            if ('*' in line):
                break
            for t in line.split(','):
                try:
                    nl.append(int(t)-1) # -1 to start nunbering from 0 
                except ValueError:
                    pass
        set_list.extend(nl)
    return set_list