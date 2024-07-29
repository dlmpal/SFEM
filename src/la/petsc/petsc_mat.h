#pragma once

#ifdef SFEM_HAS_PETSC

#include "../../common/config.h"
#include <petscmat.h>
#include <vector>

namespace sfem::la::petsc
{
    /// @brief Thin wrapper around the PETSc Mat
    class PetscMat
    {
    public:
        /// @brief Create a PetscMat
        /// @param diag_nnz Number of non-zeros for rows on the diagonal
        /// @param off_diag_nnz Number of non-zeros for rows on the off-diagonal
        PetscMat(const std::vector<int> &diag_nnz, const std::vector<int> &off_diag_nnz);

        /// @brief Create a PetscMat from an existing PETSc Mat
        /// @param A Existing PETSc Mat
        /// @param inc_ref_count Whether to increase the ref count for A
        PetscMat(Mat mat, bool inc_ref_count);

        // Copy constructor (deleted)
        PetscMat(const PetscMat &) = delete;

        // Copy assignment operator (deleted)
        PetscMat &operator=(const PetscMat &) = delete;

        /// @brief Move constructor
        PetscMat(PetscMat &&);

        /// @brief Move assignment
        PetscMat &operator=(PetscMat &&);

        /// @brief Destructor
        ~PetscMat();

        /// @brief Get the local size
        int size_local() const;

        /// @brief Get the global size
        int size_global() const;

        /// @brief Get the underlying PETSc Mat
        Mat mat() const;

        /// @brief Reset the preallocation
        /// @note Call before re-assembling
        void reset();

        /// @brief Add values to the matrix
        /// @param idxs Indices
        /// @param values Values
        void add_values(const std::vector<int> &idxs, const std::vector<Scalar> &values);

        /// @brief Assemble the matrix
        void assemble();

    private:
        /// @brief Underlying PETSc Mat
        Mat mat_;
    };
}

#endif // SFEM_HAS_PETSC