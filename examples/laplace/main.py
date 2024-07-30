"""
Solve Laplace's equation in 2D, modelling potential flow over a quarter cylinder.
"""

import sys
import petsc4py
import pysfem

# Init PETSc
petsc4py.init(sys.argv)

# Read mesh
mesh = pysfem.io.read_mesh("mesh")

# Material proprties
prop = pysfem.fe.constitutive.ThermoMechanicalProprties()
prop.kappa = 1.0

# Constitutive
thick = 1.0
type = pysfem.fe.constitutive.ThermoElasticPlaneConstitutive.Type.plane_stress
constitutive = pysfem.fe.constitutive.ThermoElasticPlaneConstitutive(
    prop, thick, type)

# Potential field
phi = pysfem.mesh.Field("Potential", 1, mesh, ["phi"])
phi.add_fixed_dof("Right", 0, 0)

# Create elements
elems = []
for cell in mesh.get_region_cells("Fluid"):
    elem = pysfem.fe.thermal.HeatConduction2D(cell, constitutive)
    elems.append(elem)
for cell in mesh.get_region_cells("Left"):
    elem = pysfem.fe.thermal.HeatFlux2D(cell, -1.0, thick)
    elems.append(elem)

# Create linear system matrix/vectors
K = pysfem.la.petsc.create_mat(mesh, 1)
F = pysfem.la.petsc.create_vec(mesh, 1)
U = pysfem.la.petsc.create_vec(mesh, 1)

# Assemble system
pysfem.fe.assemble_matrix(
    elems, phi, pysfem.fe.FEMatrixType.stiffness, K)
pysfem.fe.assemble_vector(elems, phi, pysfem.fe.FEVectorType.load, F)

# Apply Dirichlet B.C. and solve
pysfem.la.petsc.apply_fixed_dof(phi.get_fixed_dof(),
                                phi.get_fixed_dof_values(),
                                K, F, U)
pysfem.la.petsc.solve(K, F, U)

# Update field values
phi.set_values(U.get_values())

# Project the von Mises stress on to the FE basis
velocity = pysfem.fe.project_function(
    elems, phi, pysfem.fe.function.FieldGradient(phi))

# Write to file
pysfem.io.write_field_values("fields/phi_0", phi, True)
pysfem.io.write_field_values("fields/U_0", velocity, True)
# pysfem.io.write_vtk("fields/sfem.vtk", mesh, [phi, velocity])
