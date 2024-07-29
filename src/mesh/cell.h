#pragma once

#include "../geo/vec3.h"
#include <map>
#include <vector>

namespace sfem::mesh
{
    /// @brief Available cell types
    enum class CellType : int
    {
        point = 0,
        line = 1,
        triangle = 2,
        quad = 3,
        tet = 4,
        hex = 5,
        prism = 6
    };

    /// @brief Reference dimension for each cell type
    /// @brief If the returned value is -1, the cell type is invalid
    inline int cell_dim(CellType type)
    {
        static const std::map<CellType, int> dim_per_type = {
            {CellType::point, 0},
            {CellType::line, 1},
            {CellType::triangle, 2},
            {CellType::quad, 2},
            {CellType::tet, 3},
            {CellType::hex, 3},
            {CellType::prism, 3}};

        if (dim_per_type.count(type) == 0)
        {
            return -1;
        }
        else
        {
            return dim_per_type.at(type);
        }
    }

    /// @brief Number of nodes for each cell type and order
    /// @brief If the returned value is -1, the cell type or the order is invalid
    inline int cell_nodes(CellType type, int order)
    {
        static const std::map<std::pair<CellType, int>, int> n_nodes_per_type = {
            {{CellType::point, 1}, 1},

            {{CellType::line, 1}, 2},
            {{CellType::line, 2}, 3},
            {{CellType::line, 3}, 4},

            {{CellType::triangle, 1}, 3},
            {{CellType::triangle, 2}, 6},
            {{CellType::triangle, 3}, 10},

            {{CellType::quad, 1}, 4},
            {{CellType::quad, 2}, 9},
            {{CellType::quad, 3}, 16},

            {{CellType::tet, 1}, 4},
            {{CellType::tet, 2}, 10},
            {{CellType::tet, 3}, 20},

            {{CellType::hex, 1}, 8},
            {{CellType::hex, 2}, 20},

            {{CellType::prism, 1}, 6}};

        if (n_nodes_per_type.count({type, order}) == 0)
        {
            return -1;
        }
        else
        {
            return n_nodes_per_type.at(std::make_pair(type, order));
        }
    }

    /// @brief Number of faces each cell type
    /// @brief If the returned value is -1, the cell type is invalid
    inline int cell_n_faces(CellType type)
    {
        static const std::map<CellType, int> n_faces_per_type = {
            {CellType::point, 1},
            {CellType::line, 2},
            {CellType::triangle, 3},
            {CellType::quad, 4},
            {CellType::tet, 4},
            {CellType::hex, 6},
            {CellType::prism, 5}};

        if (n_faces_per_type.count(type) == 0)
        {
            return -1;
        }
        else
        {
            return n_faces_per_type.at(type);
        }
    }

    /// @brief Get the corresponding nodes for a cell face
    /// @note Use f_idx = -1 to get the nodes required for the cell outward normal vector (not available for 3D cells)
    /// @note If the returned std::array starts with -1, either the cell type or the face index is invalid
    inline std::array<int, 4> cell_n_nodes_face(CellType type, int f_idx)
    {
        static const std::map<std::pair<CellType, int>, std::array<int, 4>> nodes_per_type_and_face = {
            {{CellType::point, 0}, {0, -1, -1, -1}},

            {{CellType::line, 0}, {0, 1, -1, -1}},
            {{CellType::line, 1}, {1, 0, -1, -1}},
            {{CellType::line, -1}, {0, 1, -1, -1}},

            {{CellType::triangle, 0}, {0, 1, -1, -1}},
            {{CellType::triangle, 1}, {1, 2, -1, -1}},
            {{CellType::triangle, 2}, {2, 0, -1, -1}},
            {{CellType::triangle, -1}, {0, 1, 0, 2}},

            {{CellType::quad, 0}, {0, 1, -1, -1}},
            {{CellType::quad, 1}, {1, 2, -1, -1}},
            {{CellType::quad, 2}, {2, 3, -1, -1}},
            {{CellType::quad, 3}, {3, 0, -1, -1}},
            {{CellType::quad, -1}, {0, 1, 0, 2}}};

        if (nodes_per_type_and_face.count({type, f_idx}) == 0)
        {
            return {-1, -1, -1, -1};
        }
        else
        {
            return nodes_per_type_and_face.at({type, f_idx});
        }
    }

    /// @brief  Mesh cell
    class Cell
    {
    public:
        /// @brief Create a Cell
        /// @param idx Global index (must be greater than or equal to zero)
        /// @param type Geometrical type, e.g. triangle, quadrilateral etc.
        /// @param order Cell order (must be greater than zero)
        /// @param region_tag Integer tag of the cell's mesh region
        Cell(int idx, CellType type, int order, int region_tag);

        /// @brief Get the cell's global index
        int idx() const;

        /// @brief Get the cell's CellType
        CellType type() const;

        /// @brief Get the cell's order/degree
        int order() const;

        /// @brief Get the cell's region tag
        int region_tag() const;

        /// @brief Get the cell's reference dimension
        int dim() const;

        /// @brief Get the number of nodes for the cell
        int n_nodes() const;

        /// @brief Get the number of faces for the cell
        int n_faces() const;

        /// @brief Get face normal vector.
        /// @note Use f_idx = -1 to get the outward normal (not available for 3D cells)
        geo::Vec3 face_normal(int f_idx, const std::vector<Scalar> &xpts) const;

    private:
        /// @brief Global index
        int idx_;

        /// @brief Type
        CellType type_;

        /// @brief Order
        int order_;

        /// @brief Integer tag of the Region to which the Cell belongs
        int region_tag_;

        /// @brief Reference dimension
        int dim_;

        /// @brief Number of nodes
        int n_nodes_;

        /// @brief Number of faces
        int n_faces_;
    };
}