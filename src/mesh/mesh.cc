#include "mesh.h"
#include "../common/logger.h"
#include "../common/error.h"

namespace sfem::mesh
{
    //=============================================================================
    Mesh::Mesh(std::vector<Cell> cells, Connectivity conn, std::vector<Scalar> xpts,
               std::vector<Region> regions, common::IndexMap cell_im, common::IndexMap node_im)
        : cells_(cells), conn_(conn), xpts_(xpts),
          regions_(regions), node_im_(node_im), cell_im_(cell_im)
    {
        // Check for possible size mismatches
        if (cells.size() != static_cast<std::size_t>(conn.n1))
        {
            error::invalid_size_error(cells.size(), conn.n1, __FILE__, __LINE__);
        }
        if (cells.size() != static_cast<std::size_t>(cell_im.n_local()))
        {
            error::invalid_size_error(cells.size(), cell_im.n_local(), __FILE__, __LINE__);
        }
        if (xpts.size() / 3 != static_cast<std::size_t>(conn.n2))
        {
            error::invalid_size_error(xpts.size() / 3, conn.n2, __FILE__, __LINE__);
        }
        if (xpts.size() / 3 != static_cast<std::size_t>(node_im.n_local()))
        {
            error::invalid_size_error(xpts.size() / 3, node_im.n_local(), __FILE__, __LINE__);
        }

        // Compute physical dimension
        dim_ = 0;
        for (auto region : regions)
        {
            if (region.dim() > dim_)
            {
                dim_ = region.dim();
            }
        }
    }
    //=============================================================================
    void Mesh::info() const
    {
        std::string message = "Mesh Info\n";
        message += "Number of Nodes (Owned) " + std::to_string(node_im_.n_owned()) + " (Ghost) " + std::to_string(node_im_.n_ghost()) + "\n";
        message += "Number of Cells (Owned) " + std::to_string(cell_im_.n_owned()) + " (Ghost) " + std::to_string(cell_im_.n_ghost()) + "\n";
        message += "Number of Regions: " + std::to_string(regions_.size()) + "\n";
        message += "Regions: \n";
        for (auto region : regions_)
        {
            message += "\t|-" + region.name() + "\n";
        }
        Logger::instance().info(message);
    }
    //=============================================================================
    int Mesh::dim() const
    {
        return dim_;
    }
    //=============================================================================
    int Mesh::n_cells_local() const
    {
        return cell_im_.n_local();
    }
    //=============================================================================
    int Mesh::n_cells_global() const
    {
        return cell_im_.n_global();
    }
    //=============================================================================
    int Mesh::n_nodes_local() const
    {
        return node_im_.n_local();
    }
    //=============================================================================
    int Mesh::n_nodes_global() const
    {
        return node_im_.n_global();
    }
    //=============================================================================
    const std::vector<Scalar> &Mesh::xpts() const
    {
        return xpts_;
    }
    //=============================================================================
    void Mesh::set_xpts(const std::vector<Scalar> &xpts)
    {
        if (xpts_.size() != xpts.size())
        {
            error::invalid_size_error(xpts_.size(), xpts.size(), __FILE__, __LINE__);
        }
        xpts_ = xpts;
    }
    //=============================================================================
    const std::vector<Cell> &Mesh::cells() const
    {
        return cells_;
    }
    //=============================================================================
    const Connectivity &Mesh::cell_node_conn() const
    {
        return conn_;
    }
    //=============================================================================
    const std::vector<Region> &Mesh::regions() const
    {
        return regions_;
    }
    //=============================================================================
    const common::IndexMap &Mesh::cell_im() const
    {
        return cell_im_;
    }
    //=============================================================================
    const common::IndexMap &Mesh::node_im() const
    {
        return node_im_;
    }
    //=============================================================================
    Region Mesh::get_region_by_name(const std::string &region_name) const
    {
        for (auto region : regions_)
        {
            if (region.name() == region_name)
            {
                return region;
            }
        }
        Logger::instance().error("Invalid region name: " + region_name + "\n", __FILE__, __LINE__);
        return Region("", -1, -1);
    }
    //=============================================================================
    std::vector<Cell> Mesh::get_region_cells(const std::string &region_name) const
    {
        auto region = get_region_by_name(region_name);
        std::vector<Cell> region_cells;
        for (auto cell : cells_)
        {
            if (cell.region_tag() == region.tag())
            {
                region_cells.push_back(cell);
            }
        }
        return region_cells;
    }
    //=============================================================================
    std::vector<int> Mesh::get_region_nodes(const std::string &region_name) const
    {
        auto region_cells = get_region_cells(region_name);
        std::vector<int> region_nodes;
        std::vector<bool> is_node_included(n_nodes_local(), false);
        for (auto cell : region_cells)
        {
            auto _cell_nodes = get_cell_nodes(cell);
            for (auto i = 0; i < cell.n_nodes(); i++)
            {
                if (is_node_included.at(_cell_nodes[i]) == false)
                {
                    is_node_included[_cell_nodes[i]] = true;
                    region_nodes.push_back(_cell_nodes[i]);
                }
            }
        }
        return region_nodes;
    }
    //=============================================================================
    std::vector<int> Mesh::get_cell_nodes(const Cell &cell) const
    {
        std::vector<int> cell_nodes(cell.n_nodes());
        int local_idx = cell_im_.global_to_local(cell.idx());
        for (int i = 0; i < cell.n_nodes(); i++)
        {
            cell_nodes[i] = conn_.idx[conn_.ptr[local_idx] + i];
        }
        return cell_nodes;
    }
    //=============================================================================
    std::vector<Scalar> Mesh::get_cell_xpts(const Cell &cell) const
    {
        std::vector<Scalar> cell_xpts(cell.n_nodes() * 3);
        auto cell_nodes = get_cell_nodes(cell);
        for (int i = 0; i < cell.n_nodes(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                cell_xpts[i * 3 + j] = xpts_[cell_nodes[i] * 3 + j];
            }
        }
        return cell_xpts;
    }
}