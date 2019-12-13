# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 22:08:27 2019

@author: Boyang
"""
from . import Parameters as param

def vtkoutput1part(jobname, nodes, elem_lists, f, a, RF, NDIM, NDOF_NODE):
    """program to write mesh, nodal solutions, reaction forces to vtk format"""
    
    # set output file name
    outputfile = jobname+'-output.vtk'
    
    # set the no. of decimals for real outputs
    nd = param.num_decimal_for_output
    
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
        # calculate total no. of elems and total size of data and form the list
        # of element number for vtk output
        nelem = 0
        nsize = 0
        elnumlist = []
        for elist in elem_lists:
            nelem += len(elist.elems)
            if elist.eltype in param.tuple_frame2d2_eltypes:
                # frame2d2 elem has 3 data: elem index and 2 node indices
                nsize += 3*len(elist.elems)
                # frame2d2 elem has eltype number 3 in vtk format (see below)
                for i in range(len(elist.elems)): elnumlist.append(3)
            elif elist.eltype in param.tuple_tri2d3_eltypes:
                # tri elem has 4 data: elem index and node indices
                nsize += 4*len(elist.elems)
                # tri elem has eltype number 5 in vtk format (see below)
                for i in range(len(elist.elems)): elnumlist.append(5)
            else:
                print('Unsupported element type in vtkoutput: '+elist.eltype)
                
        # write line for total count of elems and integer numbers to be written
        output.write('CELLS '+str(nelem)+' '+str(nsize)+'\n')
        
        # write each element data: total no. of nodes and its nodal connec for
        # each line
        for elist in elem_lists:
            for elem in elist.elems:
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
        for elnum in elnumlist:
            output.write(str(elnum)+'\n')
            
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
        


        
def vtkoutput1part_igpoints(jobname, elem_lists, NST):
    """program to write integration point data such as its position, 
    displacement vector, strain and stress to vtk format"""
    
    # set output file name
    outputfile = jobname+'-output-igpoints.vtk'
    
    # set the no. of decimals for real outputs
    nd = param.num_decimal_for_output
    
    with open(outputfile,'w') as output:
        
        # write header
        output.write('# vtk DataFile Version 3.1\n')
        output.write('for FEM Programme igpoints output\n')
        output.write('ASCII\n')
        output.write('DATASET POLYDATA\n')
        
        # write igpoints
        
        # calculate the total number of igpoints and form all the lists for
        # output
        nsize = 0
        xlist = []
        ulist = []
        strainlist = []
        stresslist = []
        for elist in elem_lists:
            if elist.eltype in param.tuple_tri2d3_eltypes:
                # tri elem has 1 ig point
                nsize += len(elist.elems)
                for elem in elist.elems:
                    xlist.append(elem.igpoints[0].x)
                    ulist.append(elem.igpoints[0].u)
                    strainlist.append(elem.igpoints[0].strain)
                    stresslist.append(elem.igpoints[0].stress)      
            else:
                print('Unsupported element type in vtkoutput: '+elist.eltype)
        
        # write the positions of the igpoints
        output.write('POINTS '+str(nsize)+' double\n')
        for x in xlist:
            for xi in x:
                output.write(' '+str(round(xi,nd)))
            if (len(x) == 2): output.write(' 0.0')
            output.write('\n')
            
        # write the displacement vectors of the igpoints
        output.write('POINT_DATA '+str(nsize)+'\n')
        output.write('VECTORS displacement double\n')
        for u in ulist:
            for uj in u:
                output.write(str(round(uj,nd))+'  ')
            if (len(u) == 2): output.write('0.0')
            output.write('\n')
            
        # write strain and stress as field data
        nfield = 2
        
        output.write('FIELD FieldData '+str(nfield)+'\n') 
        # write applied forces and moments
        output.write('Strain '+str(NST)+' '+str(nsize)+' double\n')
        for strain in strainlist:
            for eps in strain[:,0]:
                output.write(str(round(eps,nd))+' ')
            output.write('\n')
        # write reaction forces and moments
        output.write('Stress '+str(NST)+' '+str(nsize)+' double\n')
        for stress in stresslist:
            for sig in stress[:,0]:
                output.write(str(round(sig,nd))+' ')
            output.write('\n')
            
        output.close()
        