#pragma once

#include "petsc_vec.h"
#include "petsc_mat.h"
#include "petsc_ksp.h"
#include "../sparsity_pattern.h"

namespace sfem::la::petsc
{
    /// @brief Create a PetscVec for a given mesh and number of variables per node
    inline PetscVec create_vec(const mesh::Mesh &mesh, int n_vars)
    {
        auto im = mesh.node_im().renumber();
        auto ghost_nodes = im.get_ghost_idxs();
        std::vector<int> ghost_dof(ghost_nodes.size() * n_vars);
        for (std::size_t i = 0; i < ghost_nodes.size(); i++)
        {
            for (int j = 0; j < n_vars; j++)
            {
                ghost_dof[i * n_vars + j] = ghost_nodes[i] * n_vars + j;
            }
        }
        return PetscVec(im.n_owned() * n_vars,
                        im.n_global() * n_vars,
                        ghost_dof);
    }

    /// @brief Create a PetscMat for a given mesh and number of variables per node
    inline PetscMat create_mat(const mesh::Mesh &mesh, int n_vars)
    {
        auto [diag_nnz, off_diag_nnz] = sparsity_pattern(mesh, n_vars);
        return PetscMat(diag_nnz, off_diag_nnz);
    }

    /// @brief For a linear system of the form Ax=b,
    /// remove rows and columns of A corresponding to fixed DoF,
    /// and add their contribution to b
    /// @note Does not affect the symmetry of A
    /// @note Calls .assemble() the A, b and x
    /// @param idxs Indices of the fixed DoF
    /// @param values Values of the fixed DoF
    /// @param A Left-hand-side (LHS) matrix
    /// @param b Right-hand-side (RHS) vector
    /// @param x Solution vector
    inline void apply_fixed_dof(const std::vector<int> &idxs,
                                const std::vector<Scalar> &values,
                                PetscMat &A,
                                PetscVec &b,
                                PetscVec &x)
    {
        A.assemble();
        x.insert_values(idxs, values);
        MatZeroRowsColumns(A.mat(), idxs.size(), idxs.data(), 1.0, x.vec(), b.vec());
        b.assemble();
        x.assemble();
    }

    /// @brief Solve the linear system Ax=b using PETSc's KSP solvers
    /// @param A Left-hand-side (LHS) matrix
    /// @param b Right-hand-side (RHS) vector
    /// @param x Solution vector
    /// @return Number of iterations performed
    inline int solve(const PetscMat &A,
                     const PetscVec &b,
                     PetscVec &x)
    {
        // Create and setup the linear solver
        la::petsc::PetscKSP solver;
        solver.set_from_options();
        solver.set_operator(A.mat());

        // Solve and return iteration number;
        int n_iter = solver.solve(b.vec(), x.vec());
        return n_iter;
    }
}