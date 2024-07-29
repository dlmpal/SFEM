#pragma once

#include "connectivity.h"
#include "../common/index_map.h"

namespace sfem::mesh
{
    /// @brief Base class for mesh partitioners
    class Partitioner
    {
    public:
        /// @brief Create a Partitioner
        /// @param n_parts Desired number of partitions
        /// @param cell_node_conn Cell-to-node connectivity of the mesh
        Partitioner(int n_parts, const Connectivity &cell_node_conn);

        /// @brief Destructor (virtual)
        virtual ~Partitioner() = 0;

        /// @brief
        std::pair<common::IndexMap, common::IndexMap> part_mesh() const;

    protected:
        /// @brief Relevant info for mesh partitions
        struct PartitionData
        {
            std::vector<int> cell_idxs;
            std::vector<int> node_idxs;
            std::vector<int> ghost_idxs;
            std::vector<int> ghost_owners;
        };

        /// @brief Compute the owning process of each cell/node.
        /// @note This function effectively partitions the mesh.
        /// Here the underlying partitioner, e.g METIS, is called.
        /// @return The list of the owning process for each cell/node
        virtual std::pair<std::vector<int>, std::vector<int>> compute_owners() const = 0;

        /// @brief Compute the ghost nodes for each process.
        /// @param node_owners List of owning process for each node.
        /// @param data List of PartitionData for each process.
        void compute_ghost_nodes(const std::vector<int> &node_owners, std::vector<PartitionData> &data) const;

        std::pair<common::IndexMap, common::IndexMap> distribute_partition_data(const std::vector<PartitionData> &data) const;

        /// @brief Number of partitions
        int n_parts_;

        /// @brief Cell-to-node connectivity of the mesh to be partitioned
        const Connectivity &conn_;
    };

    /// @brief Constructs a Partitioner given the type, e.g METIS
    Partitioner *create_partitioner(const std::string &type, int n_parts, const Connectivity &conn);

#ifdef SFEM_HAS_METIS
    /// @brief Partition the Mesh using METIS
    class METISPartitioner : public Partitioner
    {
    public:
        METISPartitioner(int n_parts, const Connectivity &conn);

    private:
        std::pair<std::vector<int>, std::vector<int>> compute_owners() const override;
    };
#endif // SFEM_HAS_METIS
}