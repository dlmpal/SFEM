#pragma once

#include "../mesh/mesh.h"
#include "../common/error.h"

namespace sfem::io
{
    /// @brief Read a Gmsh .msh2 file into sfem format
    /// @note This function is not designed for parallel execution
    mesh::Mesh read_gmsh(const std::string &path);
}

namespace sfem::io::gmsh
{
    /// @brief Get the sfem cell type and order from the corresponding Gmsh type
    inline std::pair<mesh::CellType, int> gmsh_type_to_native(int gmsh_id, int gmsh_type)
    {
        static const std::map<int, std::pair<mesh::CellType, int>> from_gmsh = {
            {15, {mesh::CellType::point, 1}},

            {1, {mesh::CellType::line, 1}},
            {8, {mesh::CellType::line, 2}},
            {26, {mesh::CellType::line, 3}},

            {2, {mesh::CellType::triangle, 1}},
            {9, {mesh::CellType::triangle, 2}},
            {21, {mesh::CellType::triangle, 3}},

            {3, {mesh::CellType::quad, 1}},
            {10, {mesh::CellType::quad, 2}},
            {36, {mesh::CellType::quad, 3}},

            {4, {mesh::CellType::tet, 1}},
            {11, {mesh::CellType::tet, 2}},
            {29, {mesh::CellType::tet, 3}},

            {5, {mesh::CellType::hex, 1}},
            {17, {mesh::CellType::hex, 2}},

            {6, {mesh::CellType::prism, 1}},
        };
        if (from_gmsh.count(gmsh_type) == 0)
        {
            error::unsupported_gmsh_type_error(gmsh_id, gmsh_type, __FILE__, __LINE__);
        }
        return from_gmsh.at(gmsh_type);
    }
}