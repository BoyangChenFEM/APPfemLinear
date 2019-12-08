# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 22:08:27 2019

@author: Boyang
"""
from .Elements import tri2d3elem

def vtkoutput1part(jobname, nodes, elems, f, a, RF, NDIM, NDOF_NODE):
    
    # set output file name
    outputfile = jobname+'-output.vtk'
    
    # set the no. of decimals for real outputs
    nd = 6
    
    with open(outputfile,'w') as output:
        
        # write header
        output.write('# vtk DataFile Version 3.1\n')
        output.write('for FEM Programme output\n')
        output.write('ASCII\n')
        output.write('DATASET UNSTRUCTURED_GRID\n')
        
        # write nodes
        output.write('POINTS '+str(len(nodes))+' double\n')
        for node in nodes:
            for x in node:
                output.write(' '+str(round(x,nd)))
            if (len(node) == 2): output.write(' 0.0')
            output.write('\n')
        
        # write element connec
        nelem = len(elems) # total no. of elems
        nsize = 0 # initialize total no. of integers to be written for elem cnc
        # calculate nsize
        for ie, elem in enumerate(elems):
            if isinstance(elem, tri2d3elem.tri2d3elem):
                nsize += 4
            else:
                print('Unsupported element type in vtkoutput for elem '+str(ie+1))
        # write line for total count of elems and integer numbers to be written
        output.write('CELLS '+str(nelem)+' '+str(nsize)+'\n')
        # write each element data: total no. of nodes and its nodal connec for
        # each line
        for elem in elems:
            output.write(str(len(elem.cnc_node))) # total no. of nodes
            for inode in elem.cnc_node: # nodal connectivity of this elem
                output.write(' '+str(inode))
            output.write('\n')
        
        # write the vtk eltype numbers for the different element types:
        output.write('CELL_TYPES'+' '+str(nelem)+'\n')
        # ---- see the below table for reference -------------
        #1	VTK_VERTEX	Vertex
        #2	VTK_POLY_VERTEX	Vertex
        #3	VTK_LINE	Edge Lagrange P1
        #5	VTK_TRIANGLE	Triangle Lagrange P1
        #8	VTK_PIXEL	Quadrilateral Lagrange P1
        #9	VTK_QUAD	Quadrilateral Lagrange P1
        #10	VTK_TETRA	Tetrahedron Lagrange P1
        #11	VTK_VOXEL	Hexahedron Lagrange P1
        #12	VTK_HEXAHEDRON	Hexahedron Lagrange P1
        #13	VTK_WEDGE	Wedge Lagrange P1
        #21	VTK_QUADRATIC_EDGE	Edge Lagrange P2
        #22	VTK_QUADRATIC_TRIANGLE	Triangle Lagrange P2
        #23	VTK_QUADRATIC_QUAD	Quadrilateral Lagrange P2
        #24	VTK_QUADRATIC_TETRA	Tetrahedron Lagrange P2
        #25	VTK_QUADRATIC_HEXAHEDRON	Hexahedron Lagrange P2
        # ----------------------------------------------------
        for ie, elem in enumerate(elems):
            if isinstance(elem, tri2d3elem.tri2d3elem):
                output.write(str(5)+'\n')
            else:
                print('Unsupported element type in vtkoutput for elem '+str(ie+1))
            
        # write nodal displacement vectors
        nnode = len(nodes)
        output.write('POINT_DATA '+str(nnode)+'\n')
        output.write('VECTORS displacement double\n')
        for i in range(nnode):
            for j in range(NDIM):
                output.write(str(round(a[i*NDOF_NODE+j,0],nd))+'  ')
            if (NDIM == 2): output.write('0.0')
            output.write('\n')
            
        # write field data:
        # - applied forces/moments f
        # - reaction forces/moments RF
        # - additional nodal DoFs if present
        
        # determine the number of field data
        if NDOF_NODE > NDIM:
            nfield = 3
        else:
            nfield = 2
        output.write('FIELD FieldData '+str(nfield)+'\n') 
        # write applied forces and moments
        output.write('AppliedForcesMoments '+str(NDOF_NODE)+' '+str(nnode)+' double\n')
        for i in range(nnode):
            for j in range(NDOF_NODE):
                output.write(str(round(f[i*NDOF_NODE+j,0],nd))+' ')
            output.write('\n')
        # write reaction forces and moments
        output.write('ReactionForcesMoments '+str(NDOF_NODE)+' '+str(nnode)+' double\n')
        for i in range(nnode):
            for j in range(NDOF_NODE):
                output.write(str(round(RF[i*NDOF_NODE+j,0],nd))+' ')
            output.write('\n')
        # write additional DoFs
        if (NDOF_NODE > NDIM):
            output.write('AdditionalDoFs '+str(NDOF_NODE-NDIM)+' '+str(nnode)+' double\n')
            for i in range(nnode):
                for j in range(NDIM, NDOF_NODE):
                    output.write(str(round(a[i*NDOF_NODE+j,0],nd))+' ')
                output.write('\n')
                
        output.close()