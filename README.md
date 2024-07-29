# 𝙎𝙁𝙀𝙈
SFEM is a C++ finite element framework, primarily developed for the solution of problems in solid mechanics and heat transfer.
It supports distributed memory parallelism via the MPI protocol.
SFEM utilizes PETSc for sparse linear algebra related computations, and METIS for mesh partitioning. 
SLEPc is also an optional dependency, which enables the computation of eigenvalues of discretized operators.
 
## Requirements
* C++ 17 and above
* PETSc (linear algebra)
* MPI (parallel execution)
* METIS (mesh partitioning)
* SLEPc (optional, eigenvalue computation)
* nanobind (optional, Python bindings)

## Installation
Start by cloning the repository to your machine. Inside the root directory for SFEM, run the install.py
script, which is responbile for configuration, building and installation. To inspect the available
arguments for install.py, run install.py -h. Some arguments, e.g. --petsc_dir are required.
It is highly recommended to set the path to your SFEM installation directory as an enviroment variable,
e.g. $SFEM_DIR.

## Example Programs
Example use of SFEM can be found under the examples/ folder.

## Non-native Mesh Formats
Currently, the only non-native formats supported by SFEM are Gmsh for importing,
and VTK for exporting meshes. Programs for doing so are found under apps/mesh_utils.