# APPfemLinear
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/zenodo/10.5281/zenodo.3242074/)

Linear FEM programme in Python, using Abaqus/Calculix input file, Python solver and Paraview postprocessing

Objective: 
personal use, for education, quick verification, unit testing, fun, whatever. 
It Supports multiple materials, element types and load types (concentrated & distributed) in the same problem.

Content:
- FEA2D-JOBNAME.py: Python job script which calls the input file (JOBNAME.inp), the program (LinearFEMProgram) and 
  writes the output files (JOBNAME-output/output-igpoints.vtk). Users define the basic dimension data, element types, 
  materials, loads, boundary conditions. They should match the element type, nsets, elsets in the input file.
  
- JOBNAME.inp: the Abaqus input file, single part is supported, nsets for boundary conditions and elsets for 
  different element/material types. It is used only to provide mesh, node set and element set information. 
  
- JOBNAME-output(-igpoints).vtk: vtk outputs of displacements and stress/strain (at igpoints), visualize using Paraview.

- LinearFEMProgram: contains all the source codes of the linear FEM program, job-independent.

- Examples: contains the script, input & outputs of some example problems

Run:
- run the FEA2D-JOBNAME.py using Spyder or something equivalent. It will generate the outputs automatically in this directory.
- Keep the FEA2D-JOBNAME.py and the corresponding input file in this directory when running.

How to create inp file using SALOME (free preprocessing software: https://www.salome-platform.org/)
1. create the geometry and the mesh following the standard salome procedure (see tutorial: https://www.youtube.com/watch?v=TbZDXt_VSTE)
2. create the groups for boundary conditions and loads
3. remove the boundary elements using "Modification-remove-elements-set filter-create-free boundary"
4. export into .unv format
5. download unv2ccx.exe from https://github.com/calculix/unv2ccx, it is a program to transform .unv file to .inp file
6. extract the unv2ccx.exe in a convenient folder, say under Desktop. Open cmd, in Windows type ".\unv2ccx.exe xxx.unv"
   the xx.inp will be automatically generated. It can then be used with APPfemLinear.
