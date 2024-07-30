// Solve Laplace's equation in 2D, modelling potential flow over a quarter cylinder.

#include "sfem.h"
#include <iostream>

using namespace sfem;

void solve_laplace(const std::string &mesh_path)
{
    auto mesh = io::read_mesh(mesh_path);

    mesh::Field phi("Potential", 1, mesh, {"phi"});
    phi.add_fixed_dof("Right", 0, 0);

    fe::constitutive::ThermoMechanicalProperties prop;
    prop.kappa = 1.0;
    Scalar thick = 1.0;
    fe::constitutive::ThermoElasticPlaneConstitutive constitutive(prop, thick);

    std::vector<std::shared_ptr<fe::FiniteElement>> elems;
    for (const auto &cell : mesh.get_region_cells("Fluid"))
    {
        auto ele = std::make_shared<fe::thermal::HeatConduction2D>(cell, constitutive);
        elems.push_back(std::move(ele));
    }
    for (const auto &cell : mesh.get_region_cells("Left"))
    {
        auto ele = std::make_shared<fe::thermal::HeatFlux2D>(cell, -1.0, 1.0);
        elems.push_back(std::move(ele));
    }

    auto K = la::petsc::create_mat(mesh, 1);
    auto U = la::petsc::create_vec(mesh, 1);
    auto F = la::petsc::create_vec(mesh, 1);

    fe::assemble_matrix(elems, phi, fe::FEMatrixType::stiffness, K);
    fe::assemble_vector(elems, phi, fe::FEVectorType::load, F);

    la::petsc::apply_fixed_dof(phi.get_fixed_dof(),
                               phi.get_fixed_dof_values(),
                               K, F, U);
    la::petsc::solve(K, F, U);

    phi.set_values(U.get_values());
    io::write_field_values("fields/phi_0", phi);

    auto velocity = fe::project_function(elems, phi, fe::function::FieldGradient(phi));
    io::write_field_values("fields/U_0", velocity);
}

int main(int argc, char **argv)
{
    initialize(&argc, &argv, "LaplaceSolver");

    std::string mesh_path = argv[1];
    solve_laplace(mesh_path);

    finalize();
    return 0;
}