#include "petsc_mat.h"
#include "../../common/logger.h"

namespace sfem::la::petsc
{
    //=============================================================================
    PetscMat::PetscMat(const std::vector<int> &diag_nnz, const std::vector<int> &off_diag_nnz)
    {
        MatCreateAIJ(SFEM_COMM_WORLD,
                     diag_nnz.size(), diag_nnz.size(),
                     PETSC_DECIDE, PETSC_DECIDE,
                     PETSC_DECIDE, diag_nnz.data(),
                     PETSC_DECIDE, off_diag_nnz.data(),
                     &mat_);

        MatSetFromOptions(mat_);
    }
    //=============================================================================
    PetscMat::PetscMat(Mat mat, bool inc_ref_count)
        : mat_(mat)
    {
        if (!mat_)
        {
            Logger::instance().error("Invalid PETSc Mat passed to PetscMat constructor", __FILE__, __LINE__);
        }

        if (inc_ref_count)
        {
            PetscObjectReference((PetscObject)mat_);
        }
    }
    //=============================================================================
    PetscMat::PetscMat(PetscMat &&other)
    {
        mat_ = other.mat_;
        other.mat_ = nullptr;
    }
    //=============================================================================
    PetscMat &PetscMat::operator=(PetscMat &&other)
    {
        if (this != &other)
        {
            if (mat_)
            {
                MatDestroy(&mat_);
            }
            mat_ = other.mat_;
            other.mat_ = nullptr;
        }

        return *this;
    }
    //=============================================================================
    PetscMat::~PetscMat()
    {
        if (mat_)
        {
            MatDestroy(&mat_);
        }
    }
    //=============================================================================
    int PetscMat::size_local() const
    {
        PetscInt n;
        MatGetLocalSize(mat_, &n, nullptr);
        return n;
    }
    //=============================================================================
    int PetscMat::size_global() const
    {
        PetscInt n;
        MatGetSize(mat_, &n, nullptr);
        return n;
    }
    //=============================================================================
    Mat PetscMat::mat() const { return mat_; }
    //=============================================================================
    void PetscMat::reset()
    {
        MatResetPreallocation(mat_);
    }
    //=============================================================================
    void PetscMat::add_values(const std::vector<int> &idxs, const std::vector<Scalar> &values)
    {
        MatSetValues(mat_,
                     idxs.size(), idxs.data(),
                     idxs.size(), idxs.data(),
                     values.data(), ADD_VALUES);
    }
    //=============================================================================
    void PetscMat::assemble()
    {
        MatAssemblyBegin(mat_, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat_, MAT_FINAL_ASSEMBLY);
    }
}