#include "petsc_ksp.h"
#include "../../common/logger.h"
#include "../../common/timer.h"

namespace sfem::la::petsc
{
    //=============================================================================
    PetscKSP::PetscKSP()
    {
        KSPCreate(SFEM_COMM_WORLD, &ksp_);
    }
    //=============================================================================
    PetscKSP::PetscKSP(KSP ksp, bool inc_ref_count)
        : ksp_(ksp)
    {
        if (!ksp_)
        {
            Logger::instance().error("Invalid PETSc KSP passed to PetscKSP constructor", __FILE__, __LINE__);
        }

        if (inc_ref_count)
        {
            PetscObjectReference((PetscObject)ksp_);
        }
    }
    //=============================================================================
    PetscKSP::PetscKSP(PetscKSP &&other)
    {
        ksp_ = other.ksp_;
        other.ksp_ = nullptr;
    }
    //=============================================================================
    PetscKSP &PetscKSP::operator=(PetscKSP &&other)
    {
        if (this != &other)
        {
            if (ksp_)
            {
                KSPDestroy(&ksp_);
            }
            ksp_ = other.ksp_;
            other.ksp_ = nullptr;
        }

        return *this;
    }
    //=============================================================================
    PetscKSP::~PetscKSP()
    {
        if (ksp_)
        {
            KSPDestroy(&ksp_);
        }
    }
    //=============================================================================
    KSP PetscKSP::ksp() const
    {
        return ksp_;
    }
    //=============================================================================
    void PetscKSP::set_from_options() const
    {
        KSPSetFromOptions(ksp_);
    }
    //=============================================================================
    void PetscKSP::set_options_prefix(const std::string &prefix) const
    {
        KSPSetOptionsPrefix(ksp_, prefix.c_str());
    }
    //=============================================================================
    void PetscKSP::set_operator(const Mat A) const
    {
        KSPSetOperators(ksp_, A, A);
    }
    //=============================================================================
    int PetscKSP::solve(const Vec b, Vec x) const
    {
        common::Timer timer("PetscKSP");

        KSPSolve(ksp_, b, x);

        // Get the number of iterations the KSP performed
        PetscInt n_iter;
        KSPGetIterationNumber(ksp_, &n_iter);

        // Check for convergence
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp_, &reason);
        if (reason < 0)
        {
            std::string message = "PetscKSP did not converge in " + std::to_string(n_iter) + " iterations\n";
            message += "Reason: " + std::to_string(reason) + "\n";
            Logger::instance().warn(message, __FILE__, __LINE__);
        }

        return static_cast<int>(n_iter);
    }
}