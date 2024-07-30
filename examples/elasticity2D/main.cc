// Solve the equations of linear elasticity in 2D.

#include "sfem.h"
#include <iostream>
using namespace sfem;

void solve_elasticity(const std::string &mesh_path, Scalar E, Scalar nu, Scalar thick)
{
    auto mesh = io::read_mesh(mesh_path);
    mesh.info();

    // Displacement field
    mesh::Field disp("U", 2, mesh, {"u", "v"});
    disp.add_fixed_dof("Fixed", 0, 0);
    disp.add_fixed_dof("Fixed", 1, 0);

    // Pressure field
    mesh::Field P("P", 1, mesh);
    P.set_all({1000});

    // Matieral proprties and constitutive law
    fe::constitutive::ThermoMechanicalProperties prop;
    prop.E = E;
    prop.nu = nu;
    prop.rho = 1e3;
    prop.alpha = 1e-5;
    fe::constitutive::ThermoElasticPlaneConstitutive constitutive(prop, thick);

    // Create elements
    std::vector<std::shared_ptr<fe::FiniteElement>> solid_elems;
    for (const auto &cell : mesh.get_region_cells("Solid"))
    {
        auto elem = std::make_unique<fe::solid::LinearElasticity2D>(cell, constitutive);
        solid_elems.push_back(std::move(elem));
    }
    std::vector<std::shared_ptr<fe::FiniteElement>> boundary_elems;
    for (const auto &cell : mesh.get_region_cells("Left"))
    {
        auto elem = std::make_unique<fe::solid::PressureLoad2D>(cell, P, thick);
        boundary_elems.push_back(std::move(elem));
    }

    // Compute structural mass
    auto mass = fe::assemble_function(solid_elems, disp, fe::function::StructuralMass2D());
    Logger::instance().info("Structural mass: " + std::to_string(mass.at(0, 0)) + "\n");

    auto K = la::petsc::create_mat(mesh, 2);
    auto F = la::petsc::create_vec(mesh, 2);
    auto U = la::petsc::create_vec(mesh, 2);

    // Assemble system: K U = F
    fe::assemble_matrix(solid_elems, disp, fe::FEMatrixType::stiffness, K);
    fe::assemble_vector(solid_elems, disp, fe::FEVectorType::load, F);
    fe::assemble_vector(boundary_elems, disp, fe::FEVectorType::load, F);

    // Remove Dirichlet B.C. and solve system
    la::petsc::apply_fixed_dof(disp.get_fixed_dof(),
                               disp.get_fixed_dof_values(),
                               K, F, U);
    la::petsc::solve(K, F, U);

    // Update field values and write to file
    disp.set_values(U.get_values());
    io::write_field_values("fields/U_0", disp);

    // Compute gradient and write to file
    auto stress = fe::project_function(solid_elems, disp, fe::function::Stress2D());
    io::write_field_values("fields/stress_0", stress);
}

int main(int argc, char **argv)
{
    initialize(&argc, &argv, "ElasticitySolver2D");

    std::string mesh_path = argv[1];
    Scalar E = 5e9;
    Scalar nu = 0.35;
    Scalar thick = 1e-3;

    solve_elasticity(mesh_path, E, nu, thick);

    finalize();
    return 0;
}