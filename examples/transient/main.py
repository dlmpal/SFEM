"""
Solve a transient heat transfer problem.
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
prop.kappa = 27
prop.cp = 420
prop.rho = 8040

# Constitutive
thick = 1.0
type = pysfem.fe.constitutive.ThermoElasticPlaneConstitutive.Type.plane_stress
constitutive = pysfem.fe.constitutive.ThermoElasticPlaneConstitutive(
    prop, thick, type)

# Temperature field
temp = pysfem.mesh.Field("Temperature", 1, mesh, ["T"])
temp.set_all([0.0])
temp.add_fixed_dof("Fixed", 0, 0.0)
temp.add_fixed_dof("Upper", 0, 100)

# Create elements
elems = []
for cell in mesh.get_region_cells("Solid"):
    elem = pysfem.fe.thermal.HeatConduction2D(cell, constitutive)
    elems.append(elem)

# Time step
dt = 0.01

# Number of time steps
Nt = 100

for i in range(Nt):
    print(f"Time: {i * dt} [s]")

    # Create linear system matrix/vectors
    M = pysfem.la.petsc.create_mat(mesh, 1)
    K = pysfem.la.petsc.create_mat(mesh, 1)
    F = pysfem.la.petsc.create_vec(mesh, 1)
    U = pysfem.la.petsc.create_vec(mesh, 1)

    # Update/set the solution vector values from the field
    U.insert_values(temp.get_local_dof(), temp.values())

    # Assemble system
    pysfem.fe.assemble_matrix(elems, temp, pysfem.fe.FEMatrixType.mass, M)
    pysfem.fe.assemble_matrix(
        elems, temp, pysfem.fe.FEMatrixType.stiffness, K)
    pysfem.fe.assemble_vector(elems, temp, pysfem.fe.FEVectorType.load, F)

    # Apply Dirichlet B.C. to system
    pysfem.la.petsc.apply_fixed_dof(temp.get_fixed_dof(),
                                    temp.get_fixed_dof_values(),
                                    K, F, U)

    # Implicit Euler timestepping
    # (M/dt + K) U_new = M/dt U_old + F
    pysfem.la.petsc.mat_scale(1/dt, M)
    pysfem.la.petsc.mat_axpy(1, M, K)
    pysfem.la.petsc.mat_mult_add(M, U, F, F)
    pysfem.la.petsc.solve(K, F, U)

    # Update field values
    temp.set_values(U.get_values())

    # Write to file
    pysfem.io.write_field_values(f"fields/T_{i}", temp, True)
    # pysfem.io.write_vtk(f"fields/sfem_{i}.vtk", mesh, [temp])
