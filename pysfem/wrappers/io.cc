#include "sfem.h"
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

using namespace sfem::io;
namespace nb = nanobind;
using namespace nb::literals;

namespace sfem_wrappers
{
    void init_io(nb::module_ &m)
    {
        // Mesh
        m.def("read_mesh", &sfem::io::read_mesh, "dir"_a, "partitioner_type"_a = "METIS");
        m.def("write_mesh", &sfem::io::write_mesh);

        // Gmsh
        m.def("read_gmsh", &sfem::io::read_gmsh);

        // VTK
        m.def("write_vtk", &sfem::io::write_vtk);

        // Field
        m.def("read_field_values", &sfem::io::read_field_values);
        m.def("write_field_values", &sfem::io::write_field_values);
    }
}