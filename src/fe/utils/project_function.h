#pragma once

#include "../finite_element.h"
#include "../functions/function.h"
#include "../../la/petsc/petsc_utils.h"
#include "../../mesh/field.h"

namespace sfem::fe
{
    inline mesh::Field project_function(const std::vector<std::shared_ptr<FiniteElement>> &elems,
                                        const mesh::Field &field,
                                        const Function &func)
    {

        // Resulting field
        mesh::Field projection("Projection",
                               func.size(),
                               field.mesh());

        auto M = la::petsc::create_mat(field.mesh(), 1); // LHS projection/mass matrix
        std::vector<la::petsc::PetscVec> F;              // RHS vector(s)
        std::vector<la::petsc::PetscVec> U;              // Solution vector(s)
        for (int i = 0; i < func.size(); i++)
        {
            F.push_back(la::petsc::create_vec(field.mesh(), 1));
            U.push_back(la::petsc::create_vec(field.mesh(), 1));
        }

        // Assemble the matrix and vector(s)
        for (const auto &elem : elems)
        {
            auto xpts = field.mesh().get_cell_xpts(elem->cell());
            auto dof = field.dof_im().local_to_global(field.mesh().get_cell_nodes(elem->cell()));
            auto u = field.get_cell_values(elem->cell());

            auto [Me, Fe] = elem->project_function(xpts, u, func);

            M.add_values(dof, Me.entries());
            for (int i = 0; i < func.size(); i++)
            {
                F[i].add_values(dof, Fe.get_col_values(i));
            }
        }

        // Solve the resulting linear system(s)
        M.assemble();
        for (int i = 0; i < func.size(); i++)
        {
            F[i].assemble();
            U[i].assemble();
            la::petsc::solve(M, F[i], U[i]);
        }

        // Set the projection field's values
        std::vector<Scalar> values(field.mesh().n_nodes_local() * func.size());
        for (int i = 0; i < func.size(); i++)
        {
            auto u = U[i].get_values();
            for (std::size_t j = 0; j < u.size(); j++)
            {
                values[j * func.size() + i] = u[j];
            }
        }
        projection.set_values(values);

        return projection;
    }
}