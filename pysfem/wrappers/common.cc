#include "sfem.h"
#include <nanobind/nanobind.h>

using namespace sfem::common;
namespace nb = nanobind;

namespace sfem_wrappers
{
    void init_common(nb::module_ &m)
    {
        // IndexMap
        nb::class_<IndexMap>(m, "IndexMap")
            .def(nb::init<int>())
            .def(nb::init<const std::vector<int> &, const std::vector<int> &, const std::vector<int> &>())
            .def("n_owned", &IndexMap::n_owned)
            .def("n_ghost", &IndexMap::n_ghost)
            .def("n_local", &IndexMap::n_local)
            .def("n_global", &IndexMap::n_global)
            .def("get_owned_idxs", &IndexMap::get_owned_idxs, nb::rv_policy::reference_internal)
            .def("get_ghost_idxs", &IndexMap::get_ghost_idxs, nb::rv_policy::reference_internal)
            .def("get_ghost_owners", &IndexMap::get_ghost_owners, nb::rv_policy::reference_internal)
            .def("local_to_global", nb::overload_cast<int>(&IndexMap::local_to_global, nb::const_))
            .def("global_to_local", nb::overload_cast<int>(&IndexMap::global_to_local, nb::const_))
            .def("renumber", &IndexMap::renumber);
    }
}