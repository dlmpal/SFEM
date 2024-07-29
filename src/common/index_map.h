#pragma once

#include "config.h"
#include <vector>
#include <unordered_map>

namespace sfem::common
{
    class IndexMap
    {
    public:
        /// @brief Create an IndexMap for serial execution
        /// @param n_owned Number of locally owned indices
        IndexMap(int n_owned);

        /// @brief Create an IndexMap
        /// @param owned_idxs Locally owned indices
        /// @param ghost_idxs Ghost indices for this process
        /// @param ghost_owners Owning process for each ghost index
        IndexMap(const std::vector<int> &owned_idxs, const std::vector<int> &ghost_idxs, const std::vector<int> &ghost_owners);

        /// @brief Get the number of indices owned by this process
        int n_owned() const;

        /// @brief Get the number of ghost indices for this process
        int n_ghost() const;

        /// @brief Get the number of local indices for this process
        /// @note Sum of owned and ghost indices
        int n_local() const;

        /// @brief Get the global number of indices
        int n_global() const;

        /// @brief Get the owned indices
        std::vector<int> get_owned_idxs() const;

        /// @brief Get the ghost indices
        std::vector<int> get_ghost_idxs() const;

        /// @brief Get the owning process for each ghost index
        std::vector<int> get_ghost_owners() const;

        /// @brief Map an index or set of indices from local to global indexing
        /// @note Raises an error if idx > n_local
        int local_to_global(int idx) const;

        /// @brief Apply local_to_global(int) to an array of indices
        std::vector<int> local_to_global(const std::vector<int> &idxs) const;

        /// @brief Map an index from global to local indexing
        /// @note Returns -1 if the index is not local
        int global_to_local(int idx) const;

        /// @brief Apply global_to_local(int) to an array of indices
        std::vector<int> global_to_local(const std::vector<int> &idxs) const;

        /// @brief
        /// @return
        IndexMap renumber() const;

    private:
        int n_owned_ = 0;
        int n_ghost_ = 0;
        std::vector<int> local_to_global_;
        std::unordered_map<int, int> global_to_local_;
        std::vector<int> ghost_owners_;
    };
}