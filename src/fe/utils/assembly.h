#pragma once

#include "../finite_element.h"
#include "../../la/petsc/petsc_mat.h"
#include "../../la/petsc/petsc_vec.h"
#include "../../mesh/field.h"
#include "../../common/timer.h"

namespace sfem::fe
{

    inline void assemble_matrix(const std::vector<std::shared_ptr<FiniteElement>> &elems,
                                const mesh::Field &field,
                                FEMatrixType type,
                                la::petsc::PetscMat &mat)
    {
        // Time the assembly
        common::Timer timer("Matrix assembly");

        auto &mesh = field.mesh();

        for (const auto &elem : elems)
        {
            // Cell data
            auto xpts = mesh.get_cell_xpts(elem->cell());
            auto dof = field.get_cell_dof(elem->cell());
            auto u = field.get_cell_values(elem->cell());

            // Integrate and add contribution
            auto elem_matrix = elem->integrate_fe_matrix(xpts, u, type);
            mat.add_values(dof, elem_matrix.entries());
        }
    }

    inline void assemble_vector(const std::vector<std::shared_ptr<FiniteElement>> &elems,
                                const mesh::Field &field,
                                FEVectorType type,
                                la::petsc::PetscVec &vec)
    {
        // Time the assembly
        common::Timer timer("Vector assembly");

        auto &mesh = field.mesh();

        for (const auto &elem : elems)
        {
            // Cell data
            auto xpts = mesh.get_cell_xpts(elem->cell());
            auto dof = field.get_cell_dof(elem->cell());
            auto u = field.get_cell_values(elem->cell());

            // Integrate and add contribution
            auto elem_vec = elem->integrate_fe_vector(xpts, u, type);
            vec.add_values(dof, elem_vec.entries());
        }
    }

    /// @brief Assemble (integrate) a function for the given elements
    /// @param elems Elements to use for integration
    /// @param field Corresponding field
    /// @param func Function to be integrated
    /// @return The integrated function value(s)
    inline la::DenseMatrix assemble_function(const std::vector<std::shared_ptr<FiniteElement>> &elems,
                                             const mesh::Field &field,
                                             const Function &func)
    {
        la::DenseMatrix value_(func.size(), 1);

        for (const auto &elem : elems)
        {
            // Cell data
            auto xpts = field.mesh().get_cell_xpts(elem->cell());
            auto dof = field.get_cell_dof(elem->cell());
            auto u = field.get_cell_values(elem->cell());

            value_ += elem->integrate_function(xpts, u, func);
        }

        if (Logger::instance().n_procs() > 1)
        {
            std::vector<Scalar> recv_buffer(func.size());
            MPI_Allreduce(value_.entries().data(),
                          recv_buffer.data(),
                          func.size(),
                          MPI_DOUBLE,
                          MPI_SUM,
                          SFEM_COMM_WORLD);

            la::DenseMatrix value(func.size(), 1, recv_buffer);
            return value;
        }
        else
        {
            return value_;
        }
    }
}