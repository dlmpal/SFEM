// Solve the equations of linear elasticity in 3D.

#include "sfem.h"

using namespace sfem;

void solve_elasticity(const std::string &mesh_path, Scalar rho, Scalar E, Scalar nu)
{
    // Read mesh
    auto mesh = io::read_mesh(mesh_path);

    // Displacement field
    mesh::Field disp("U", 3, mesh, {"U", "V", "W"});
    disp.add_fixed_dof("Fixed", 0, 0);
    disp.add_fixed_dof("Fixed", 1, 0);
    disp.add_fixed_dof("Fixed", 2, 0);

    // Material proprties
    fe::constitutive::ThermoMechanicalProperties prop;
    prop.E = E;
    prop.nu = nu;
    prop.rho = rho;

    // Constitutive
    fe::constitutive::ThermoElasticSolidConstitutive te_const(prop);

    // Create elements
    std::vector<std::shared_ptr<fe::FiniteElement>> solid_elems;
    for (const auto &cell : mesh.get_region_cells("Solid"))
    {
        auto ele = std::make_shared<fe::solid::LinearElasticity3D>(cell, te_const);
        ele->add_inertial_load({0, -9.81, 0});
        solid_elems.push_back(std::move(ele));
    }

    // Compute structural mass
    auto mass = fe::assemble_function(solid_elems, disp, fe::function::StructuralMass3D());
    Logger::instance().info("Structural mass: " + std::to_string(mass.at(0, 0)) + "\n");

    // Create linear system matrix/vectors
    auto K = la::petsc::create_mat(mesh, 3);
    auto F = la::petsc::create_vec(mesh, 3);
    auto U = la::petsc::create_vec(mesh, 3);

    // Assemble system
    fe::assemble_matrix(solid_elems, disp, fe::FEMatrixType::stiffness, K);
    fe::assemble_vector(solid_elems, disp, fe::FEVectorType::load, F);

    // Apply Dirichlet B.C. and solve
    la::petsc::apply_fixed_dof(disp.get_fixed_dof(),
                               disp.get_fixed_dof_values(),
                               K, F, U);
    la::petsc::solve(K, F, U);

    // Update field values
    disp.set_values(U.get_values());

    // Project the von Mises stress on to the FE basis
    auto stress_vm = fe::project_function(solid_elems, disp, fe::function::VonMisesStress3D());

    // Write to file
    io::write_field_values("fields/U_0", disp);
    io::write_field_values("fields/stress_0", stress_vm);
}

int main(int argc, char **argv)
{
    initialize(&argc, &argv, "ElasticitySolver3D");

    std::string mesh_path = argv[1];
    Scalar rho = 1e-3;
    Scalar E = 1e5;
    Scalar nu = 0.2;

    solve_elasticity(mesh_path, rho, E, nu);

    finalize();
    return 0;
}