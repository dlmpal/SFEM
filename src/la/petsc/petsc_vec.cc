#include "petsc_vec.h"
#include "../../common/logger.h"

namespace sfem::la::petsc
{
    //=============================================================================
    PetscVec::PetscVec(int n_local, int n_global, const std::vector<int> &ghosts)
    {
        VecCreateGhost(SFEM_COMM_WORLD, n_local, n_global, ghosts.size(),
                       ghosts.data(), &vec_);
        set_all(0.0);
    }
    //=============================================================================
    PetscVec::PetscVec(Vec vec, bool inc_ref_count) : vec_(vec)
    {
        if (!vec_)
        {
            Logger::instance().error("Invalid PETSc Vec passed to PetscVec constructor", __FILE__, __LINE__);
        }

        if (inc_ref_count)
        {
            PetscObjectReference((PetscObject)vec_);
        }
    }
    //=============================================================================
    PetscVec::PetscVec(PetscVec &&other)
    {
        vec_ = other.vec_;
        other.vec_ = nullptr;
    }
    //=============================================================================
    PetscVec &PetscVec::operator=(PetscVec &&other)
    {
        if (this != &other)
        {
            if (vec_)
            {
                VecDestroy(&vec_);
            }
            vec_ = other.vec_;
            other.vec_ = nullptr;
        }

        return *this;
    }
    //=============================================================================
    PetscVec::~PetscVec()
    {
        if (vec_)
        {
            VecDestroy(&vec_);
        }
    }
    //=============================================================================
    int PetscVec::size_local() const
    {
        PetscInt n;
        VecGetLocalSize(vec_, &n);
        return n;
    }
    //=============================================================================
    int PetscVec::size_global() const
    {
        PetscInt n;
        VecGetSize(vec_, &n);
        return n;
    }
    //=============================================================================
    Vec PetscVec::vec() const { return vec_; }
    //=============================================================================
    PetscVec PetscVec::copy() const
    {
        Vec vec;
        VecDuplicate(vec_, &vec);
        VecCopy(vec_, vec);
        PetscVec copy(vec, false);
        return copy;
    }
    //=============================================================================
    void PetscVec::set_all(Scalar value) { VecSet(vec_, value); }
    //=============================================================================
    void PetscVec::add_values(const std::vector<int> &idxs, const std::vector<Scalar> &values)
    {
        VecSetValues(vec_, idxs.size(), idxs.data(), values.data(), ADD_VALUES);
    }
    //=============================================================================
    void PetscVec::insert_values(const std::vector<int> &idxs, const std::vector<Scalar> &values)
    {
        VecSetValues(vec_, idxs.size(), idxs.data(), values.data(), INSERT_VALUES);
    }
    //=============================================================================
    void PetscVec::assemble()
    {
        VecAssemblyBegin(vec_);
        VecAssemblyEnd(vec_);
    }
    //=============================================================================
    std::vector<Scalar> PetscVec::get_values() const
    {
        Vec x_local;
        PetscInt n;
        Scalar *_values;
        std::vector<Scalar> values;

        VecGhostUpdateBegin(vec_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(vec_, INSERT_VALUES, SCATTER_FORWARD);

        VecGhostGetLocalForm(vec_, &x_local);
        VecGetLocalSize(x_local, &n);
        VecGetArray(x_local, &_values);

        values.resize(n);
        for (int i = 0; i < n; i++)
        {
            values[i] = _values[i];
        }

        VecRestoreArray(x_local, &_values);
        VecGhostRestoreLocalForm(vec_, &x_local);

        return values;
    }
}