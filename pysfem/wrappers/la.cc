#include "sfem.h"
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>

using namespace sfem::la;
namespace nb = nanobind;

namespace sfem_wrappers
{
    void init_petsc(nb::module_ &m);

    void init_la(nb::module_ &m)
    {
        // DenseMatrix
        nb::class_<DenseMatrix>(m, "DenseMatrix")
            .def(nb::init<int, int>())
            .def(nb::init<int, int, const std::vector<Scalar> &>())
            .def("n_rows", &DenseMatrix::n_rows)
            .def("n_cols", &DenseMatrix::n_cols)
            .def("entries", &DenseMatrix::entries, nb::rv_policy::reference_internal)
            .def("set_all", &DenseMatrix::set_all)
            .def("at", &DenseMatrix::at)
            .def("insert", &DenseMatrix::insert)
            .def("add", &DenseMatrix::add)
            .def("T", &DenseMatrix::T)
            .def(nb::self * nb::self)
            .def(nb::self * Scalar())
            .def(nb::self + nb::self)
            .def(nb::self += nb::self)
            .def(nb::self + Scalar())
            .def(nb::self += Scalar())
            .def(nb::self - nb::self)
            .def(nb::self -= nb::self)
            .def(nb::self - Scalar())
            .def(nb::self -= Scalar());

        // PETSc
        nb::module_ petsc = m.def_submodule("petsc", "PETSc");
        init_petsc(petsc);
    }

    void init_petsc(nb::module_ &m)
    {
        using namespace petsc;

        // PetscVec
        nb::class_<PetscVec>(m, "PetscVec")
            .def(nb::init<int, int, const std::vector<int> &>())
            .def("size_local", &PetscVec::size_local)
            .def("size_global", &PetscVec::size_global)
            .def("copy", &PetscVec::copy)
            .def("set_all", &PetscVec::set_all)
            .def("add_values", &PetscVec::add_values)
            .def("insert_values", &PetscVec::insert_values)
            .def("assemble", &PetscVec::assemble)
            .def("get_values", &PetscVec::get_values);

        // PetscMat
        nb::class_<PetscMat>(m, "PetscMat")
            .def(nb::init<const std::vector<int> &,
                          const std::vector<int> &>())
            .def("size_local", &PetscMat::size_local)
            .def("size_global", &PetscMat::size_global)
            .def("reset", &PetscMat::reset)
            .def("local_values", &PetscMat::add_values)
            .def("assemble", &PetscMat::assemble);

        // PETSc utils
        m.def("create_vec", &create_vec);
        m.def("create_mat", &create_mat);
        m.def("vec_scale", &vec_scale);
        m.def("mat_mult_add", &mat_mult_add);
        m.def("mat_scale", &mat_scale);
        m.def("mat_axpy", &mat_axpy);
        m.def("apply_fixed_dof", &apply_fixed_dof);
        m.def("solve", &solve);

        // // PetscKSP
        // nb::class_<PetscKSP>(m, "PetscKSP")
        //     .def(nb::init<>())
        //     .def("set_from_options", &PetscKSP::SetFromOptions)
        //     .def("set_options_prefix", &PetscKSP::SetOptionsPrefix)
        //     .def("set_operators", &PetscKSP::SetOperators)
        //     .def("solve", &PetscKSP::Solve);

        // // SlepcEPS
        // nb::class_<SlepcEPS>(m, "SlepcEPS")
        //     .def(nb::init<>())
        //     .def("set_from_options", &SlepcEPS::SetFromOptions)
        //     .def("set_num_eigenpairs", &SlepcEPS::SetNumEigenpairs)
        //     //.def("set_operators", nb::overload_cast<const PetscMat &, const PetscMat &>(&SlepcEPS::SetOperators))
        //     //.def("set_operators", nb::overload_cast<const PetscMat &>(&SlepcEPS::SetOperators))
        //     .def("solve", &PetscKSP::Solve)
        //     .def("get_num_eigenpairs", &SlepcEPS::GetNumEigenPairs)
        //     .def("get_eigenvalue", &SlepcEPS::GetEigenvalue)
        //     .def("get_eigepair", &SlepcEPS::GetEigenpair);
    }
}