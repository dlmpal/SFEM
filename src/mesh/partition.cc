#include "partition.h"
#include "../common/logger.h"
#include "../common/error.h"
#include <mpi.h>

#ifdef SFEM_HAS_METIS
#include <metis.h>
#endif

namespace sfem::mesh
{
    //=============================================================================
    Partitioner::Partitioner(int n_parts, const Connectivity &conn)
        : n_parts_(n_parts), conn_(conn)
    {
    }
    //=============================================================================
    Partitioner::~Partitioner()
    {
    }
    //=============================================================================
    std::pair<common::IndexMap, common::IndexMap> Partitioner::part_mesh() const
    {
        // The root process partitions the mesh
        std::vector<PartitionData> data_per_proc(n_parts_);
        if (Logger::instance().proc_rank() == SFEM_ROOT)
        {

            // Compute the cell and node owners for each process
            auto [cell_owners, node_owners] = compute_owners();

            // Count owned cells per process
            std::vector<int> n_cells_per_proc(n_parts_, 0);
            for (int i : cell_owners)
            {
                n_cells_per_proc[i]++;
            }

            // Assign owned cells to each process
            for (int i = 0; i < n_parts_; i++)
            {
                data_per_proc[i].cell_idxs.reserve(n_cells_per_proc[i]);
            }
            for (std::size_t i = 0; i < cell_owners.size(); i++)
            {
                data_per_proc[cell_owners[i]].cell_idxs.push_back(i);
            }

            // Count owned nodes per process
            std::vector<int> n_nodes_per_proc(n_parts_, 0);
            for (int i : node_owners)
            {
                n_nodes_per_proc[i]++;
            }

            // Assign owned nodes to each process
            for (int i = 0; i < n_parts_; i++)
            {
                data_per_proc[i].node_idxs.reserve(n_nodes_per_proc[i]);
            }
            for (std::size_t i = 0; i < node_owners.size(); i++)
            {
                data_per_proc[node_owners[i]].node_idxs.push_back(i);
            }

            // Compute the ghost nodes for each process
            compute_ghost_nodes(node_owners, data_per_proc);
        }

        // Data is distributed to all processes by root
        return distribute_partition_data(data_per_proc);
    }
    //=============================================================================
    void Partitioner::compute_ghost_nodes(const std::vector<int> &node_owners, std::vector<PartitionData> &data_per_proc) const
    {
        // Count ghost nodes per process
        std::vector<int> n_ghosts_per_proc(n_parts_, 0);
        for (int i = 0; i < n_parts_; i++)
        {
            std::unordered_map<int, short> is_node_included;
            for (std::size_t j = 0; j < data_per_proc[i].cell_idxs.size(); j++)
            {
                int cell_idx = data_per_proc[i].cell_idxs[j];
                for (int k = 0; k < conn_.cnt[cell_idx]; k++)
                {
                    int node_idx = conn_.idx[conn_.ptr[cell_idx] + k];
                    int owner = node_owners[node_idx];
                    if (is_node_included.count(node_idx) == 0 && owner != i)
                    {
                        is_node_included[node_idx] = 1;
                        n_ghosts_per_proc[i]++;
                    }
                }
            }
        }

        // Assign ghost nodes to each process
        for (int i = 0; i < n_parts_; i++)
        {
            data_per_proc[i].ghost_idxs.reserve(n_ghosts_per_proc[i]);
            data_per_proc[i].ghost_owners.reserve(n_ghosts_per_proc[i]);
        }
        for (int i = 0; i < n_parts_; i++)
        {
            std::unordered_map<int, short> is_node_included;
            for (std::size_t j = 0; j < data_per_proc[i].cell_idxs.size(); j++)
            {
                int cell_idx = data_per_proc[i].cell_idxs[j];
                for (int k = 0; k < conn_.cnt[cell_idx]; k++)
                {
                    int node_idx = conn_.idx[conn_.ptr[cell_idx] + k];
                    int owner = node_owners[node_idx];
                    if (is_node_included.count(node_idx) == 0 && owner != i)
                    {
                        is_node_included[node_idx] = 1;
                        data_per_proc[i].ghost_idxs.push_back(node_idx);
                        data_per_proc[i].ghost_owners.push_back(owner);
                    }
                }
            }
        }
    }
    //=============================================================================
    std::pair<common::IndexMap, common::IndexMap> Partitioner::distribute_partition_data(const std::vector<PartitionData> &data_per_proc) const
    {
        PartitionData data;

        // Concatenate all partition-related vectors
        // Cells
        std::vector<int> cell_idxs_per_proc;
        std::vector<int> cell_counts(n_parts_, 0);
        std::vector<int> cell_displs(n_parts_, 0);
        // Nodes
        std::vector<int> node_idxs_per_proc;
        std::vector<int> node_counts(n_parts_, 0);
        std::vector<int> node_displs(n_parts_, 0);
        // Ghosts
        std::vector<int> ghost_idxs_per_proc;
        std::vector<int> ghost_owners_per_proc;
        std::vector<int> ghost_counts(n_parts_, 0);
        std::vector<int> ghost_displs(n_parts_, 0);
        if (Logger::instance().proc_rank() == SFEM_ROOT)
        {
            for (int i = 0; i < n_parts_; i++)
            {
                cell_counts[i] = data_per_proc[i].cell_idxs.size();
                node_counts[i] = data_per_proc[i].node_idxs.size();
                ghost_counts[i] = data_per_proc[i].ghost_idxs.size();
                if (i > 0)
                {
                    cell_displs[i] = cell_displs[i - 1] + data_per_proc[i - 1].cell_idxs.size();
                    node_displs[i] = node_displs[i - 1] + data_per_proc[i - 1].node_idxs.size();
                    ghost_displs[i] = ghost_displs[i - 1] + data_per_proc[i - 1].ghost_idxs.size();
                }
                cell_idxs_per_proc.reserve(cell_idxs_per_proc.size() + data_per_proc[i].cell_idxs.size());
                node_idxs_per_proc.reserve(node_idxs_per_proc.size() + data_per_proc[i].node_idxs.size());
                ghost_idxs_per_proc.reserve(ghost_idxs_per_proc.size() + data_per_proc[i].ghost_idxs.size());
                ghost_owners_per_proc.reserve(ghost_idxs_per_proc.size() + data_per_proc[i].ghost_owners.size());
            }
            for (int i = 0; i < n_parts_; i++)
            {
                for (std::size_t j = 0; j < data_per_proc[i].cell_idxs.size(); j++)
                {
                    cell_idxs_per_proc.push_back(data_per_proc[i].cell_idxs[j]);
                }
                for (std::size_t j = 0; j < data_per_proc[i].node_idxs.size(); j++)
                {
                    node_idxs_per_proc.push_back(data_per_proc[i].node_idxs[j]);
                }
                for (std::size_t j = 0; j < data_per_proc[i].ghost_idxs.size(); j++)
                {
                    ghost_idxs_per_proc.push_back(data_per_proc[i].ghost_idxs[j]);
                    ghost_owners_per_proc.push_back(data_per_proc[i].ghost_owners[j]);
                }
            }
        }

        int n_cells;
        MPI_Scatter(cell_counts.data(), 1, MPI_INT, &n_cells, 1, MPI_INT, SFEM_ROOT, SFEM_COMM_WORLD);
        data.cell_idxs.resize(n_cells);
        MPI_Scatterv(cell_idxs_per_proc.data(), cell_counts.data(), cell_displs.data(), MPI_INT, data.cell_idxs.data(), data.cell_idxs.size(), MPI_INT, SFEM_ROOT, SFEM_COMM_WORLD);

        int n_nodes;
        MPI_Scatter(node_counts.data(), 1, MPI_INT, &n_nodes, 1, MPI_INT, SFEM_ROOT, SFEM_COMM_WORLD);
        data.node_idxs.resize(n_nodes);
        MPI_Scatterv(node_idxs_per_proc.data(), node_counts.data(), node_displs.data(), MPI_INT, data.node_idxs.data(), data.node_idxs.size(), MPI_INT, SFEM_ROOT, SFEM_COMM_WORLD);

        int n_ghosts;
        MPI_Scatter(ghost_counts.data(), 1, MPI_INT, &n_ghosts, 1, MPI_INT, SFEM_ROOT, SFEM_COMM_WORLD);
        data.ghost_idxs.resize(n_ghosts);
        data.ghost_owners.resize(n_ghosts);
        MPI_Scatterv(ghost_idxs_per_proc.data(), ghost_counts.data(), ghost_displs.data(), MPI_INT, data.ghost_idxs.data(), data.ghost_idxs.size(), MPI_INT, SFEM_ROOT, SFEM_COMM_WORLD);
        MPI_Scatterv(ghost_owners_per_proc.data(), ghost_counts.data(), ghost_displs.data(), MPI_INT, data.ghost_owners.data(), data.ghost_owners.size(), MPI_INT, SFEM_ROOT, SFEM_COMM_WORLD);

        return std::make_pair(common::IndexMap(data.cell_idxs, {}, {}), common::IndexMap(data.node_idxs, data.ghost_idxs, data.ghost_owners));
    }
    //=============================================================================
    Partitioner *create_partitioner(const std::string &type, int n_parts, const Connectivity &conn)
    {
        Partitioner *partitioner = nullptr;

        if (type == "METIS")
        {
#ifdef SFEM_HAS_METIS
            partitioner = new METISPartitioner(n_parts, conn);
#else
            Logger::GetInstance().Error("SFEM was not compiled with METIS. Add SFEM_USE_METIS=ON and re-compile the library.\n", __FILE__, __LINE__);
#endif // SFEM_USE_METIS
        }

        if (partitioner == nullptr)
        {
            Logger::instance().error("Partitioner type: " + type + " is not available.\n", __FILE__, __LINE__);
        }

        return partitioner;
    }
    //=============================================================================
#ifdef SFEM_HAS_METIS
    METISPartitioner::METISPartitioner(int n_parts, const Connectivity &conn) : Partitioner(n_parts, conn)
    {
    }
    //=============================================================================
    std::pair<std::vector<int>, std::vector<int>> METISPartitioner::compute_owners() const
    {
        // METIS-required arrays
        // The eptr array contains the first and last node of each cell
        // The eind array contains the nodes corresponding to each cell
        // For example for element i eind[eptr[i] + j] is the node index of the j-th node for cell i
        // The cell_owner and node_owner arrays contain map the cells/nodes to their owning process
        int n_cells = conn_.n1, n_nodes = conn_.n2;
        int size = static_cast<int>(conn_.idx.size());
        std::vector<idx_t> eptr(n_cells + 1);
        std::vector<idx_t> eind(size);
        std::vector<idx_t> cell_owners(n_cells);
        std::vector<idx_t> node_owners(n_nodes);

        for (int i = 0; i < n_cells; i++)
        {
            eptr[i] = conn_.ptr[i];

            for (int j = 0; j < conn_.cnt[i]; j++)
            {
                eind[eptr[i] + j] = conn_.idx[conn_.ptr[i] + j];
            }
        }
        eptr[n_cells] = size;

        // Execute METIS partitioning routine
        int metis_obj_val;
        int metis_return_val = METIS_PartMeshNodal((idx_t *)&n_cells,
                                                   (idx_t *)&n_nodes,
                                                   eptr.data(),
                                                   eind.data(),
                                                   nullptr,
                                                   nullptr,
                                                   (idx_t *)&n_parts_,
                                                   nullptr,
                                                   nullptr,
                                                   (idx_t *)&metis_obj_val,
                                                   cell_owners.data(),
                                                   node_owners.data());

        // Check return value for errors
        if (metis_return_val != 1)
        {
            std::string message = "METIS_PartMeshNodal() returned with:" + std::to_string(metis_return_val) + ". Exiting\n";
            Logger::instance().error(message, __FILE__, __LINE__);
        }

        return std::make_pair(cell_owners, node_owners);
    }
#endif // SFEM_HAS_METIS
}