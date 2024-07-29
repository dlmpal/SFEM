#pragma once

#include "connectivity.h"
#include "cell.h"
#include "region.h"
#include "../common/index_map.h"

namespace sfem::mesh
{
    class Mesh
    {
    public:
        /// @brief Create a Mesh
        /// @param cells Cells
        /// @param conn Cell-to-node connectivity
        /// @param xpts Nodal coordinates
        /// @param regions Regions
        /// @param cell_im Cell IndexMap
        /// @param node_im Node IndexMap
        Mesh(std::vector<Cell> cells, Connectivity conn, std::vector<Scalar> xpts,
             std::vector<Region> regions, common::IndexMap cell_im, common::IndexMap node_im);

        /// @brief Copy constructor
        Mesh(const Mesh &) = default;

        // Copy assignment (deleted)
        Mesh &operator=(const Mesh &) = delete;

        /// @brief Move constructor
        Mesh(Mesh &&) = default;

        /// @brief Move assignment
        Mesh &operator=(Mesh &&mesh) = default;

        /// @brief Destructor
        ~Mesh() = default;

        /// @brief Print information about the Mesh
        void info() const;

        /// @brief Get the physical dimension
        int dim() const;

        /// @brief Get the local number of cells
        int n_cells_local() const;

        /// @brief Get the global number of cells
        int n_cells_global() const;

        /// @brief Get the local number of nodes
        int n_nodes_local() const;

        /// @brief Get the global number of nodes
        int n_nodes_global() const;

        /// @brief Get a reference to the Cell vector
        const std::vector<Cell> &cells() const;

        /// @brief Get a reference to the cell-to-node connectivity
        const Connectivity &cell_node_conn() const;

        /// @brief Get a reference to the nodal positions vector
        const std::vector<Scalar> &xpts() const;

        /// @brief Set the mesh nodal positions
        void set_xpts(const std::vector<Scalar> &xpts);

        /// @brief Get a reference to the Region vector
        const std::vector<Region> &regions() const;

        /// @brief Get a reference to the Cell IndexMap
        const common::IndexMap &cell_im() const;

        /// @brief Get a reference to the node IndexMap
        const common::IndexMap &node_im() const;

        /// @brief Get a Region by its name
        Region get_region_by_name(const std::string &region_name) const;

        /// @brief Get the Cells for a Region, given its name
        std::vector<Cell> get_region_cells(const std::string &region_name) const;

        /// @brief Get the nodes for a Region, given its name
        /// @note Nodes are returned in local indexing
        std::vector<int> get_region_nodes(const std::string &region_name) const;

        /// @brief Get the nodes of a cell
        /// @note Nodes are returned in local indexing
        std::vector<int> get_cell_nodes(const Cell &cell) const;

        /// @brief Get the nodal positions of a cell
        std::vector<Scalar> get_cell_xpts(const Cell &cell) const;

    private:
        /// @brief Vector containing the local Cells
        std::vector<Cell> cells_;

        /// @brief Cell-to-node connectivity
        /// @note Nodes are stored in local indexing
        Connectivity conn_;

        /// @brief Points
        std::vector<Scalar> xpts_;

        /// @brief Mesh regions
        std::vector<Region> regions_;

        /// @brief Node IndexMap
        common::IndexMap node_im_;

        /// @brief Cell IndexMap
        common::IndexMap cell_im_;

        /// @brief Physical dimension
        int dim_;
    };
}