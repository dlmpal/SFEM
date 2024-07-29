#pragma once

#include "../mesh/field.h"
#include "../common/error.h"

namespace sfem::io
{
    /// @brief Creates .vtk file for the given Mesh and Fields
    /// @note This function will create a file for each process that runs it.
    /// Therefore for each partition a separate file is created.
    void write_vtk(const std::string &path, const mesh::Mesh &mesh, const std::vector<mesh::Field> &fields);
}

namespace sfem::io::vtk
{
    /// @brief Get the corresponding VTK type for a sfem cell type
    inline int cell_type_to_vtk(const mesh::Cell &cell)
    {
        static const std::map<std::pair<mesh::CellType, int>, int> to_vtk = {
            {{mesh::CellType::point, 1}, 1},

            {{mesh::CellType::line, 1}, 3},
            {{mesh::CellType::line, 2}, 21},
            {{mesh::CellType::line, 3}, 68},

            {{mesh::CellType::triangle, 1}, 5},
            {{mesh::CellType::triangle, 2}, 22},
            {{mesh::CellType::triangle, 3}, 69},

            {{mesh::CellType::quad, 1}, 9},
            {{mesh::CellType::quad, 2}, 23},
            {{mesh::CellType::quad, 3}, 70},

            {{mesh::CellType::tet, 1}, 10},
            {{mesh::CellType::tet, 2}, 24},
            {{mesh::CellType::tet, 3}, 71},

            {{mesh::CellType::hex, 1}, 12},
            {{mesh::CellType::hex, 2}, 25},
            {{mesh::CellType::hex, 3}, 72},

            {{mesh::CellType::prism, 1}, 13}

        };

        if (to_vtk.count({cell.type(), cell.order()}) == 0)
        {
            error::invalid_cell_error(cell.idx(), static_cast<int>(cell.type()), cell.order(), __FILE__, __LINE__);
        }

        return to_vtk.at({cell.type(), cell.order()});
    }

    /// @brief Get the corresponding node ordering for VTK cells.
    /// @note For most sfem cell types no re-ordering is needed.
    inline void cell_node_ordering_to_vtk(const mesh::Cell &cell, int nodes[])
    {
        if (cell.type() == mesh::CellType::tet and cell.order() == 2)
        {
            int temp = nodes[9];
            nodes[9] = nodes[8];
            nodes[8] = temp;
        }
    }
}