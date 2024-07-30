import sys
import petsc4py
import pysfem

# Init PETSc
petsc4py.init(sys.argv)

# Read mesh
mesh = pysfem.io.read_mesh("mesh/mesh_quad2")

# Material proprties
prop = pysfem.fe.constitutive.ThermoMechanicalProprties()
prop.E = 5e9
prop.nu = 0.35

# Constitutive
thick = 1e-3
type = pysfem.fe.constitutive.ThermoElasticPlaneConstitutive.Type.plane_stress
constitutive = pysfem.fe.constitutive.ThermoElasticPlaneConstitutive(
    prop, thick, type)

# Displacement field
disp = pysfem.mesh.Field("Displacement", 2, mesh, ["Ux", "Uy"])
disp.add_fixed_dof("Fixed", 0, 0)
disp.add_fixed_dof("Fixed", 1, 0)

# Pressure field
pressure = pysfem.mesh.Field("Pressure", 1, mesh, ["P"])
pressure.set_all([1000])

# Create elements
solid_elems = []
for cell in mesh.get_region_cells("Solid"):
    elem = pysfem.fe.solid.LinearElasticity2D(cell, constitutive)
    solid_elems.append(elem)

boundary_elems = []
for cell in mesh.get_region_cells("Left"):
    elem = pysfem.fe.solid.PressureLoad2D(cell, pressure, thick)
    boundary_elems.append(elem)

# Create linear system matrix/vectors
K = pysfem.la.petsc.create_mat(mesh, 2)
F = pysfem.la.petsc.create_vec(mesh, 2)
U = pysfem.la.petsc.create_vec(mesh, 2)

# Assemble system
pysfem.fe.assemble_matrix(
    solid_elems, disp, pysfem.fe.FEMatrixType.stiffness, K)
pysfem.fe.assemble_vector(solid_elems, disp, pysfem.fe.FEVectorType.load, F)
pysfem.fe.assemble_vector(boundary_elems, disp, pysfem.fe.FEVectorType.load, F)

# Apply Dirichlet B.C. and solve
pysfem.la.petsc.apply_fixed_dof(disp.get_fixed_dof(),
                                disp.get_fixed_dof_values(),
                                K, F, U)
pysfem.la.petsc.solve(K, F, U)

# Update field values
disp.set_values(U.get_values())

# Project the von Mises stress on to the FE basis
von_mises = pysfem.fe.project_function(
    solid_elems, disp, pysfem.fe.function.Stress2D())

# Write to file
# pysfem.io.write_field_values("fields/U_0", disp, True)
# pysfem.io.write_field_values("fields/stress_0", von_mises, True)
pysfem.io.write_vtk("fields/sfem.vtk", mesh, [disp, pressure, von_mises])
