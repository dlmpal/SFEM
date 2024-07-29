import sys
import petsc4py
import pysfem

# Init PETSc
petsc4py.init(sys.argv)

# Read mesh
mesh = pysfem.io.read_mesh("mesh/mesh_tet2")

# Material proprties
prop = pysfem.fe.constitutive.ThermoMechanicalProprties()
prop.E = 1e5
prop.nu = 0.2
prop.rho = 1e-3

# Constitutive
constitutive = pysfem.fe.constitutive.ThermoElasticSolidConstitutive(prop)

# Displacement field
disp = pysfem.mesh.Field("Displacement", 3, mesh, ["Ux", "Uy", "Uz"])
disp.add_fixed_dof("Fixed", 0, 0)
disp.add_fixed_dof("Fixed", 1, 0)
disp.add_fixed_dof("Fixed", 2, 0)

# Create elements
elems = []
for cell in mesh.get_region_cells("Solid"):
    elem = pysfem.fe.solid.LinearElasticity3D(cell, constitutive)
    elem.add_inertial_load([0, -9.81, 0])
    elems.append(elem)

# Create linear system matrix/vectors
K = pysfem.la.petsc.create_mat(mesh, 3)
F = pysfem.la.petsc.create_vec(mesh, 3)
U = pysfem.la.petsc.create_vec(mesh, 3)

# Assemble system
pysfem.fe.assemble_matrix(elems, disp, pysfem.fe.FEMatrixType.stiffness, K)
pysfem.fe.assemble_vector(elems, disp, pysfem.fe.FEVectorType.load, F)

# Apply Dirichlet B.C. and solve
pysfem.la.petsc.apply_fixed_dof(disp.get_fixed_dof(),
                                disp.get_fixed_dof_values(),
                                K, F, U)
pysfem.la.petsc.solve(K, F, U)

# Update field values
disp.set_values(U.get_values())

# Project the von Mises stress on to the FE basis
von_mises = pysfem.fe.project_function(
    elems, disp, pysfem.fe.function.VonMisesStress3D())

# Write to file
pysfem.io.write_field_values("fields/U_0", disp, True)
pysfem.io.write_field_values("fields/stress_0", von_mises, True)
