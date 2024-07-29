#include "sfem.h"
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>

using namespace sfem::geo;
namespace nb = nanobind;

namespace sfem_wrappers
{
    void init_geo(nb::module_ &m)
    {
        // Vec3
        nb::class_<Vec3>(m, "Vec3")
            .def(nb::init<>())
            .def(nb::init<Scalar, Scalar, Scalar>())
            .def(nb::init<Scalar, Scalar, Scalar, Scalar, Scalar, Scalar>())
            .def(nb::self + nb::self)
            .def(nb::self - nb::self)
            .def(nb::self * nb::self)
            .def(nb::self * Scalar())
            .def("mag", &Vec3::mag)
            .def("cross", &Vec3::cross_prod)
            .def("normalize", &Vec3::normalize)
            .def("unit_tan", &Vec3::unit_tan)
            .def("unit_norm", &Vec3::unit_norm)
            .def_rw("x", &Vec3::x);
    }
}