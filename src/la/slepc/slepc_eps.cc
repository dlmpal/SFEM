#include "slepc_eps.h"
#include "../../common/logger.h"
#include "../../common/timer.h"

#ifdef SFEM_HAS_SLEPC

namespace sfem::la::slepc
{
    //=============================================================================
    SlepcEPS::SlepcEPS()
    {
        EPSCreate(SFEM_COMM_WORLD, &eps_);
    }
    //=============================================================================
    SlepcEPS::SlepcEPS(EPS eps, bool inc_ref_count)
        : eps_(eps)
    {
        if (!eps_)
        {
            Logger::instance().error("Invalid SLEPc EPS passed to SlepcEPS constructor", __FILE__, __LINE__);
        }

        if (inc_ref_count)
        {
            PetscObjectReference((PetscObject)eps_);
        }
    }
    //=============================================================================
    SlepcEPS::SlepcEPS(SlepcEPS &&other)
    {
        eps_ = other.eps_;
        other.eps_ = nullptr;
    }
    //=============================================================================
    SlepcEPS &SlepcEPS::operator=(SlepcEPS &&other)
    {
        if (this != &other)
        {
            if (eps_)
            {
                EPSDestroy(&eps_);
            }
            eps_ = other.eps_;
            other.eps_ = nullptr;
        }

        return *this;
    }
    //=============================================================================
    SlepcEPS::~SlepcEPS()
    {
        if (eps_)
        {
            EPSDestroy(&eps_);
        }
    }
    //=============================================================================
    EPS SlepcEPS::eps() const
    {
        return eps_;
    }
    //=============================================================================
    void SlepcEPS::set_from_options() const
    {
        EPSSetFromOptions(eps_);
    }
    //=============================================================================
    void SlepcEPS::set_operators(const Mat A, const Mat B) const
    {
        EPSSetOperators(eps_, A, B);
        EPSSetProblemType(eps_, EPS_GHEP);
    }
    //=============================================================================
    void SlepcEPS::set_operators(const Mat A) const
    {
        EPSSetOperators(eps_, A, nullptr);
        EPSSetProblemType(eps_, EPS_HEP);
    }
    //=============================================================================
    int SlepcEPS::solve(int n) const
    {
        common::Timer timer("SlepcEPS");

        // TODO Check n
        EPSSetDimensions(eps_, (PetscInt)n, PETSC_DECIDE, PETSC_DECIDE);

        EPSSolve(eps_);

        // Get the number of iterations the EPS performed
        PetscInt n_iter;
        EPSGetIterationNumber(eps_, &n_iter);

        // Check for convergence
        EPSConvergedReason reason;
        EPSGetConvergedReason(eps_, &reason);
        if (reason < 0)
        {
            std::string message = "SlepcEPS did not converge in " + std::to_string(n_iter) + " iterations\n";
            message += "Reason: " + std::to_string(reason) + "\n";
            Logger::instance().warn(message, __FILE__, __LINE__);
        }

        return static_cast<int>(n_iter);
    }
    //=============================================================================
    int SlepcEPS::get_n_converged() const
    {
        PetscInt n_pairs;
        EPSGetConverged(eps_, &n_pairs);
        return static_cast<int>(n_pairs);
    }
    //=============================================================================
    std::tuple<Scalar, Scalar> SlepcEPS::get_eigenvalue(int idx) const
    {
        // Check if idx is greater than the number of converged pairs
        if (idx >= get_n_converged())
        {
            std::string message = "Eigenvalue " + std::to_string(idx) + " is not available\n";
            Logger::instance().error(message, __FILE__, __LINE__);
        }

        // Get the eigenvalue
        PetscScalar eig_real, eig_imag;
        EPSGetEigenpair(eps_, idx, &eig_real, &eig_imag, nullptr, nullptr);

        return std::make_tuple(eig_real, eig_imag);
    }
    //=============================================================================
    std::tuple<Scalar, Scalar, petsc::PetscVec, petsc::PetscVec> SlepcEPS::get_eigenpair(int idx) const
    {
        // Check if idx is greater than the number of converged pairs
        if (idx >= get_n_converged())
        {
            std::string message = "Eigenpair " + std::to_string(idx) + " is not available\n";
            Logger::instance().error(message, __FILE__, __LINE__);
        }

        // Get operator local and global size
        Mat mat;
        PetscInt n_local, n_global;
        EPSGetOperators(eps_, &mat, nullptr);
        MatGetSize(mat, &n_global, &n_local);

        // Get the eigenvalue and its corresponding eigenvector
        PetscScalar eig_real, eig_imag;
        petsc::PetscVec V_real(n_local, n_global, {});
        petsc::PetscVec V_imag(n_local, n_global, {});
        EPSGetEigenpair(eps_, idx, &eig_real, &eig_imag, V_real.vec(), V_imag.vec());

        return std::make_tuple(eig_real, eig_imag, std::move(V_real), std::move(V_imag));
    }
}

#endif // SFEM_HAS_SLEPC