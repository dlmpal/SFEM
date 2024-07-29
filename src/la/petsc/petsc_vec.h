#pragma once

#ifdef SFEM_HAS_PETSC

#include "../../common/config.h"
#include <petscvec.h>
#include <vector>

namespace sfem::la::petsc
{
    /// @brief Thin wrapper around the PETSc Vec
    class PetscVec
    {
    public:
        /// @brief Create a PetscVec
        /// @param n_local Local size
        /// @param n_global Global size
        /// @param ghosts Ghost indeces
        PetscVec(int n_local, int n_global, const std::vector<int> &ghosts);

        /// @brief Create a PetscVec from an existing PETSc Vec
        /// @param x Existing PETSc Vec
        /// @param inc_ref_count Whether to increase the reference count for x
        PetscVec(Vec x, bool inc_ref_count);

        // Copy constructor (deleted)
        PetscVec(const PetscVec &) = delete;

        // Copy assignment operator (deleted)
        PetscVec &operator=(const PetscVec &) = delete;

        /// @brief Move constructor
        PetscVec(PetscVec &&);

        /// @brief Move assignment
        PetscVec &operator=(PetscVec &&);

        /// @brief Destructor
        ~PetscVec();

        /// @brief Get the local size
        int size_local() const;

        /// @brief Get the global size
        int size_global() const;

        /// @brief Get the underlying Vec
        Vec vec() const;

        /// @brief Copy the vector
        PetscVec copy() const;

        /// @brief Set all vector values
        void set_all(Scalar value);

        /// @brief Add values into the vector
        /// @param idxs Indices
        /// @param values Values
        void add_values(const std::vector<int> &idxs, const std::vector<Scalar> &values);

        /// @brief Insert values into the vector, overriding existing ones
        /// @param idxs Indices
        /// @param values Values
        void insert_values(const std::vector<int> &idxs, const std::vector<Scalar> &values);

        /// @brief Assemble the vector
        void assemble();

        /// @brief Get the (local) values
        /// @note The values for ghost indices are also included
        std::vector<Scalar> get_values() const;

    private:
        /// @brief Underlying PETSc Vec
        Vec vec_;
    };
}

#endif // SFEM_HAS_PETSC