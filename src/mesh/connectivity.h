#pragma once

#include "../common/config.h"
#include <vector>

namespace sfem::mesh
{
    /// @brief Connectivity of n1 entities to n2 entities,
    /// e.g. how the n1 entities are connected to the n2 entities
    /// @note n1 represents the size of ptr and cnt
    /// @note n2 represents the number of unique entities in idx
    /// @note Accessing the related/connected entities of an entity i
    /// is done as follows: idx[ptr[i] + j] for j = 0, cnt[i]
    struct Connectivity
    {
        int n1 = 0;
        int n2 = 0;
        std::vector<int> ptr;
        std::vector<int> cnt;
        std::vector<int> idx;
    };

    /// @brief Invert the given connectivity, e.g get the connectivity
    // from n2 entities to n1 entities
    /// @note This is can used, for example, to get the node-to-cell connectivity
    /// given the the cell-to-node connectivity
    /// @todo Check validity
    Connectivity invert_conn(const Connectivity &conn);

    /// @brief Get the node-to-node connectivity, given the cell-to-node connectivity
    /// @todo Change to the more general "n2-to-n2" connectivity
    /// @todo Check validity
    Connectivity compute_node_to_node_conn(const Connectivity &cell_node_conn);
}