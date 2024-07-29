#pragma once

#ifdef SFEM_HAS_SLEPC

#include "../petsc/petsc_vec.h"
#include <slepc.h>
#include <tuple>

namespace sfem::la::slepc
{
    /// @brief Thin wrapper around the SLEPc eigenproblem solvers
    class SlepcEPS
    {
    public:
        /// @brief Create a SLEPcEPS object
        SlepcEPS();

        /// @brief Create a PetscMat from an existing SLEPc EPS
        /// @param ksp Existing SLEPc EPS
        /// @param inc_ref_count Whether to increase the ref count for eps
        SlepcEPS(EPS eps, bool inc_ref_count);

        // Copy constructor (deleted)
        SlepcEPS(const SlepcEPS &) = delete;

        // Copy assignment operator (deleted)
        SlepcEPS &operator=(SlepcEPS &) = delete;

        /// @brief Move constructor
        SlepcEPS(SlepcEPS &&);

        /// @brief Move assignment
        SlepcEPS &operator=(SlepcEPS &&);

        /// @brief Destructor
        ~SlepcEPS();

        /// @brief Get the underlying SLEPc EPS
        EPS eps() const;

        /// @brief Set the EPS options from the options database
        void set_from_options() const;

        /// @brief Set the operators for the generalized eigenproblem
        void set_operators(const Mat A, const Mat B) const;

        /// @brief Set the operator for the standard eigenproblem
        void set_operators(const Mat A) const;

        /// @brief Solve the eigenvalue problem using the EPS
        /// @param n Number of eigenpairs to compute
        int solve(int n) const;

        /// @brief Get the number of converged eigenpairs
        int get_n_converged() const;

        /// @brief Get the i-th eigenvalue
        std::tuple<Scalar, Scalar> get_eigenvalue(int idx) const;

        /// @brief Get the i-th eigenpair
        std::tuple<Scalar, Scalar, petsc::PetscVec, petsc::PetscVec> get_eigenpair(int idx) const;

    private:
        EPS eps_;
    };
}

#endif // SFEM_HAS_SLEPC