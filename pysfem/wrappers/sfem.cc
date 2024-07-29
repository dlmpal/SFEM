#include "sfem.h"
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>

namespace nb = nanobind;

namespace sfem_wrappers
{
    void init_common(nb::module_ &m);
    void init_geo(nb::module_ &m);
    void init_mesh(nb::module_ &m);
    void init_io(nb::module_ &m);
    void init_la(nb::module_ &m);
    void init_fe(nb::module_ &m);
};

NB_MODULE(pysfem, m)
{
    // Doc
    m.doc() = "SFEM Python interface";

    // Bind vectors
    nb::bind_vector<std::vector<Scalar>>(m, "ScalarVector");
    nb::bind_vector<std::vector<int>>(m, "IntVector");
    nb::bind_vector<std::vector<std::string>>(m, "StringVector");
    nb::bind_vector<std::vector<sfem::mesh::Cell>>(m, "CellVector");
    nb::bind_vector<std::vector<sfem::mesh::Region>>(m, "RegionVector");
    nb::bind_vector<std::vector<std::shared_ptr<sfem::fe::FiniteElement>>>(m, "ElementVector");

    // Common
    nb::module_ common = m.def_submodule("common", "Commonly used utilities");
    sfem_wrappers::init_common(common);

    // Geo
    nb::module_ geo = m.def_submodule("geo", "Computational geometry functionality");
    sfem_wrappers::init_geo(geo);

    // Mesh
    nb::module_ mesh = m.def_submodule("mesh", "Mesh-related functionality and utilities");
    sfem_wrappers::init_mesh(mesh);

    // LA
    nb::module_ la = m.def_submodule("la", "Linear Algebra");
    sfem_wrappers::init_la(la);

    // I/O
    nb::module_ io = m.def_submodule("io", "I/O utilities");
    sfem_wrappers::init_io(io);

    // FE
    nb::module_ fe = m.def_submodule("fe", "FEM");
    sfem_wrappers::init_fe(fe);
}