# APPfemLinear
Linear FEM programme in Python, using Abaqus/Salome input file, Python solver and Paraview postprocessing

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
