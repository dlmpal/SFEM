#include "sfem.h"
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

using namespace sfem::mesh;
using namespace sfem::common;
namespace nb = nanobind;

namespace sfem_wrappers
{
    void init_mesh(nb::module_ &m)
    {
        // CellType
        nb::enum_<CellType>(m, "CellType")
            .value("point", CellType::point)
            .value("line", CellType::line)
            .value("triangle", CellType::triangle)
            .value("quad", CellType::quad)
            .value("tet", CellType::tet)
            .value("hex", CellType::hex)
            .value("prism", CellType::prism);

        // Cell
        nb::class_<Cell>(m, "Cell")
            .def(nb::init<int, CellType, int, int>())
            .def("idx", &Cell::idx)
            .def("type", &Cell::type)
            .def("order", &Cell::order)
            .def("region_tag", &Cell::region_tag)
            .def("dim", &Cell::dim)
            .def("n_nodes", &Cell::n_nodes)
            .def("n_faces", &Cell::n_faces)
            .def("face_normal", &Cell::face_normal);

        // Region
        nb::class_<Region>(m, "Region")
            .def(nb::init<const std::string &, int, int>())
            .def("name", &Region::name)
            .def("dim", &Region::dim)
            .def("tag", &Region::tag);

        // Connectivity
        nb::class_<Connectivity>(m, "Connectivity")
            .def(nb::init<>())
            .def_rw("n1", &Connectivity::n1)
            .def_rw("n2", &Connectivity::n2)
            .def_rw("ptr", &Connectivity::ptr)
            .def_rw("cnt", &Connectivity::cnt)
            .def_rw("idx", &Connectivity::idx);

        // Connectivity-related utility functions
        m.def("invert_connectivity", &invert_conn);
        m.def("compute_cell_to_node_connectivity", &compute_node_to_node_conn);

        // Mesh
        nb::class_<Mesh>(m, "Mesh")
            .def(nb::init<std::vector<Cell>, Connectivity,
                          std::vector<double>, std::vector<Region>,
                          IndexMap, IndexMap>())

            .def("info", &Mesh::info)
            .def("dim", &Mesh::dim)
            .def("n_cells_local", &Mesh::n_cells_local)
            .def("n_cells_global", &Mesh::n_cells_global)
            .def("n_nodes_local", &Mesh::n_nodes_local)
            .def("n_nodes_global", &Mesh::n_nodes_global)

            .def("cells", &Mesh::cells, nb::rv_policy::reference_internal)
            .def("cell_node_conn", &Mesh::cell_node_conn, nb::rv_policy::reference_internal)
            .def("xpts", &Mesh::xpts, nb::rv_policy::reference_internal)
            .def("regions", &Mesh::regions, nb::rv_policy::reference_internal)
            .def("cell_im", &Mesh::cell_im, nb::rv_policy::reference_internal)
            .def("node_im", &Mesh::node_im, nb::rv_policy::reference_internal)

            .def("get_region_by_name", &Mesh::get_region_by_name)
            .def("get_region_cells", &Mesh::get_region_cells)
            .def("get_region_nodes", &Mesh::get_region_nodes)
            .def("get_cell_nodes", &Mesh::get_cell_nodes)
            .def("get_cell_xpts", &Mesh::get_cell_xpts);

        // Field
        nb::class_<Field>(m, "Field")
            .def(nb::init<const std::string &,
                          int,
                          Mesh &,
                          const std::vector<std::string> &>())
            .def("name", &Field::name)
            .def("n_vars", &Field::n_vars)
            .def("mesh", &Field::mesh, nb::rv_policy::reference_internal)
            .def("comp_names", &Field::comp_names)
            .def("dof_im", &Field::dof_im)
            .def("n_dof_owned", &Field::n_dof_owned)
            .def("n_dof_ghost", &Field::n_dof_ghost)
            .def("n_dof_local", &Field::n_dof_local)
            .def("n_dof_global", &Field::n_dof_global)
            .def("map_node_dof", &Field::map_node_dof)
            .def("get_owned_dof", &Field::get_owned_dof)
            .def("get_ghost_dof", &Field::get_ghost_dof)
            .def("get_cell_dof", &Field::get_cell_dof)
            .def("get_cell_values", &Field::get_cell_values)
            .def("add_fixed_dof", &Field::add_fixed_dof)
            .def("get_fixed_dof", &Field::get_fixed_dof)
            .def("get_fixed_dof_values", &Field::get_fixed_dof_values)
            .def("clear_fixed_dof", &Field::clear_fixed_dof)
            .def("set_all", &Field::set_all)
            .def("set_values", &Field::set_values)
            .def("values", &Field::values, nb::rv_policy::reference_internal);
    }
}