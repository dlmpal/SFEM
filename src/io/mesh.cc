#include "mesh.h"
#include "../common/error.h"

namespace sfem::io
{
    //=============================================================================
    std::pair<std::vector<mesh::Cell>, mesh::Connectivity>
    read_cells(const std::string &path, bool distributed, const common::IndexMap &cell_im, const common::IndexMap &node_im)
    {
        if (distributed == false && Logger::instance().proc_rank() != SFEM_ROOT)
        {
            return std::make_pair(std::vector<mesh::Cell>(), mesh::Connectivity());
        }

        std::ifstream file(path);
        if (!file.is_open())
        {
            error::invalid_filename_error(path, __FILE__, __LINE__);
        }

        // Global sizes
        int n_cells_global, n_nodes_global, conn_size_global;
        file >> n_cells_global;
        file >> n_nodes_global;
        file >> conn_size_global;

        // Local sizes
        int n_cells_local, n_nodes_local, conn_size_local;
        if (distributed == false)
        {
            n_cells_local = n_cells_global;
            n_nodes_local = n_nodes_global;
            conn_size_local = conn_size_global;
        }
        else
        {
            // Check that global sizes match
            if (cell_im.n_global() != n_cells_global)
            {
                error::invalid_size_error(cell_im.n_global(), n_cells_global, __FILE__, __LINE__);
            }
            if (node_im.n_global() != n_nodes_global)
            {
                error::invalid_size_error(node_im.n_global(), n_nodes_global, __FILE__, __LINE__);
            }

            // Local sizes known from the index maps
            n_cells_local = cell_im.n_local();
            n_nodes_local = node_im.n_local();

            // Since the local connectivity size is not known,
            // we have to loop over the cells to compute it
            conn_size_local = 0;
            for (int i = 0; i < n_cells_global; i++)
            {
                int idx_global, type, order, region_tag, node;

                file >> idx_global;
                file >> type;
                file >> order;
                file >> region_tag;

                mesh::Cell cell(idx_global, static_cast<mesh::CellType>(type), order, region_tag);

                for (int j = 0; j < cell.n_nodes(); j++)
                {
                    file >> node;
                }

                if (cell_im.global_to_local(cell.idx()) >= 0)
                {
                    conn_size_local += cell.n_nodes();
                }
            }

            // Close the file, re-open and skip the first line
            file.close();
            file.open(path);
            file >> n_cells_global;
            file >> n_nodes_global;
            file >> conn_size_global;
        }

        // Read cells and cell-to-node connectivity
        std::vector<mesh::Cell> cells;
        cells.reserve(n_cells_local);
        mesh::Connectivity conn;
        conn.n1 = n_cells_local;
        conn.n2 = n_nodes_local;
        conn.ptr.resize(n_cells_local);
        conn.cnt.resize(n_cells_local);
        conn.idx.resize(conn_size_local);

        int ptr = 0, idx_local = 0;
        std::vector<int> cell_nodes;
        for (int i = 0; i < n_cells_global; i++)
        {
            int idx_global, type, order, region_tag;

            file >> idx_global;
            file >> type;
            file >> order;
            file >> region_tag;

            mesh::Cell cell(idx_global, static_cast<mesh::CellType>(type), order, region_tag);

            if (static_cast<int>(cell_nodes.size()) != cell.n_nodes())
            {
                cell_nodes.resize(cell.n_nodes());
            }

            for (int j = 0; j < cell.n_nodes(); j++)
            {
                file >> cell_nodes[j];
            }

            // For distributed meshes, add only if cell is local to this process
            if (distributed == false || cell_im.global_to_local(idx_global) >= 0)
            {
                cells.push_back(cell);
                conn.ptr[idx_local] = ptr;
                conn.cnt[idx_local] = cell.n_nodes();
                for (int j = 0; j < cell.n_nodes(); j++)
                {
                    conn.idx[ptr + j] = cell_nodes[j];
                }
                ptr += cell.n_nodes();
                idx_local++;
            }
        }

        // For distributed meshes, map the nodes to local indexing
        if (distributed)
        {
            for (std::size_t i = 0; i < conn.idx.size(); i++)
            {
                conn.idx[i] = node_im.global_to_local(conn.idx[i]);
            }
        }

        return std::make_pair(cells, conn);
    }
    //=============================================================================
    std::vector<Scalar> read_xpts(const std::string &path, bool distributed, const common::IndexMap &node_im)
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            error::invalid_filename_error(path, __FILE__, __LINE__);
        }

        // Number of nodes (global)
        int n_nodes_global;
        file >> n_nodes_global;

        // Number of nodes (local)
        int n_nodes_local;
        if (distributed == false)
        {
            n_nodes_local = n_nodes_global;
        }
        else
        {
            if (node_im.n_global() != n_nodes_global)
            {
                error::invalid_size_error(node_im.n_global(), n_nodes_global, __FILE__, __LINE__);
            }

            // Local size known from the index map
            n_nodes_local = node_im.n_local();
        }

        // Read nodal coordinates
        std::vector<Scalar> xpts(n_nodes_local * 3);
        Scalar xpt[3];
        for (int i = 0; i < n_nodes_global; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                file >> xpt[j];
            }

            // For distributed meshes only add only if the node is local to this process
            if (distributed == false || node_im.global_to_local(i) >= 0)
            {
                int local_idx = distributed == false ? i : node_im.global_to_local(i);
                for (int j = 0; j < 3; j++)
                {
                    xpts[local_idx * 3 + j] = xpt[j];
                }
            }
        }

        return xpts;
    }
    //=============================================================================
    std::vector<mesh::Region> read_regions(const std::string &path)
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            error::invalid_filename_error(path, __FILE__, __LINE__);
        }

        // Read the number of regions
        int n_regions;
        file >> n_regions;

        // Read regions
        std::vector<mesh::Region> regions;
        for (auto i = 0; i < n_regions; i++)
        {
            std::string name;
            int dim, tag;
            file >> name;
            file >> dim;
            file >> tag;
            regions.push_back(mesh::Region(name, dim, tag));
        }

        return regions;
    }
    //=============================================================================
    mesh::Mesh read_mesh(const std::string &dir, const std::string &partitioner_type)
    {
        int n_procs = Logger::instance().n_procs();
        bool distributed = false;
        common::IndexMap cell_im(0);
        common::IndexMap node_im(0);

        // For distributed meshes, first create the partitions
        if (n_procs > 1)
        {
            distributed = true;
            auto [_, conn] = read_cells(dir + "/cells", false, common::IndexMap(0), common::IndexMap(0));
            auto partitioner = mesh::create_partitioner(partitioner_type, n_procs, conn);
            std::tie(cell_im, node_im) = partitioner->part_mesh();
            delete partitioner;
        }

        auto [cells, conn] = read_cells(dir + "/cells", distributed, cell_im, node_im);
        if (n_procs == 1)
        {
            cell_im = common::IndexMap(conn.n1);
            node_im = common::IndexMap(conn.n2);
        }
        auto xpts = read_xpts(dir + "/xpts", distributed, node_im);
        auto regions = read_regions(dir + "/regions");

        return mesh::Mesh(cells, conn, xpts, regions, cell_im, node_im);
    }
    //=============================================================================
    void write_cells(const std::string &path, const mesh::Mesh &mesh)
    {
        std::ofstream file(path + "/cells");
        if (!file.is_open())
        {
            error::invalid_filename_error(path + "/cells", __FILE__, __LINE__);
        }

        const auto &cells = mesh.cells();
        const auto &conn = mesh.cell_node_conn();

        file << conn.n1 << "\n";
        file << conn.n2 << "\n";
        file << conn.idx.size() << "\n";
        for (const auto &cell : cells)
        {
            file << cell.idx() << " " << static_cast<int>(cell.type()) << " " << cell.order() << " " << cell.region_tag() << " ";
            auto cell_nodes = mesh.get_cell_nodes(cell);
            for (int j = 0; j < cell.n_nodes(); j++)
            {
                file << cell_nodes[j] << " ";
            }
            file << "\n";
        }
    }
    //=============================================================================
    void write_xts(const std::string &path, const mesh::Mesh &mesh)
    {
        std::ofstream file(path + "/xpts");
        if (!file.is_open())
        {
            error::invalid_filename_error(path + "/xpts", __FILE__, __LINE__);
        }
        file << mesh.n_nodes_local() << "\n";
        const auto &xpts = mesh.xpts();
        for (int i = 0; i < mesh.n_nodes_local(); i++)
        {
            file << xpts[i * 3 + 0] << " ";
            file << xpts[i * 3 + 1] << " ";
            file << xpts[i * 3 + 2] << "\n";
        }
    }
    //=============================================================================
    void write_regions(const std::string &path, const mesh::Mesh &mesh)
    {
        std::ofstream file(path + "/regions");
        if (!file.is_open())
        {
            error::invalid_filename_error(path + "/regions", __FILE__, __LINE__);
        }

        const auto &regions = mesh.regions();

        file << regions.size() << "\n";
        for (const auto &region : regions)
        {
            file << region.name() << " " << region.dim() << " " << region.tag() << "\n";
        }
    }
    //=============================================================================
    void write_mesh(const std::string &dir, const mesh::Mesh &mesh)
    {
        write_xts(dir, mesh);
        write_cells(dir, mesh);
        write_regions(dir, mesh);
    }
}