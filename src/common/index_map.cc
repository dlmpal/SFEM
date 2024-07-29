#include "index_map.h"
#include "error.h"
#include <mpi.h>

namespace sfem::common
{
    //=============================================================================
    IndexMap::IndexMap(int n_owned)
        : n_owned_(n_owned), n_ghost_(0), ghost_owners_({})
    {
        local_to_global_.resize(n_owned);
        for (int i = 0; i < n_owned; i++)
        {
            local_to_global_[i] = i;
            global_to_local_[i] = i;
        }
    }
    //=============================================================================
    IndexMap::IndexMap(const std::vector<int> &owned_idxs, const std::vector<int> &ghost_idxs, const std::vector<int> &ghost_owners)
        : n_owned_(owned_idxs.size()), n_ghost_(ghost_idxs.size()), ghost_owners_(ghost_owners)
    {
        // Check ghost sizes match
        if (ghost_idxs.size() != ghost_owners.size())
        {
            error::invalid_size_error(ghost_idxs.size(), ghost_owners.size(), __FILE__, __LINE__);
        }

        // Create the local-to-global mapping
        local_to_global_.reserve(n_owned_ + n_ghost_);
        for (int i = 0; i < n_owned_; i++)
        {
            local_to_global_.push_back(owned_idxs[i]);
        }
        for (int i = 0; i < n_ghost_; i++)
        {
            local_to_global_.push_back(ghost_idxs[i]);
        }

        // Create the global-to-local mapping
        for (int i = 0; i < n_owned_ + n_ghost_; i++)
        {
            global_to_local_[local_to_global_.at(i)] = i;
        }
    }
    //=============================================================================
    int IndexMap::n_owned() const
    {
        return n_owned_;
    }
    //=============================================================================
    int IndexMap::n_ghost() const
    {
        return n_ghost_;
    }
    //=============================================================================
    int IndexMap::n_local() const
    {
        return n_owned_ + n_ghost_;
    }
    //=============================================================================
    int IndexMap::n_global() const
    {
        if (Logger::instance().n_procs() == 1)
        {
            return n_owned_;
        }
        else
        {
            int n_global;
            MPI_Allreduce(&n_owned_, &n_global, 1, MPI_INT, MPI_SUM, SFEM_COMM_WORLD);
            return n_global;
        }
    }
    //=============================================================================
    std::vector<int> IndexMap::get_owned_idxs() const
    {
        std::vector<int> owned_idxs(n_owned_);
        for (int i = 0; i < n_owned_; i++)
        {
            owned_idxs[i] = local_to_global_[i];
        }
        return owned_idxs;
    }
    //=============================================================================
    std::vector<int> IndexMap::get_ghost_idxs() const
    {
        std::vector<int> ghosts_idxs(n_ghost_);
        for (int i = 0; i < n_ghost_; i++)
        {
            ghosts_idxs[i] = local_to_global_[i + n_owned_];
        }
        return ghosts_idxs;
    }
    //=============================================================================
    std::vector<int> IndexMap::get_ghost_owners() const
    {
        return ghost_owners_;
    }
    //=============================================================================
    int IndexMap::local_to_global(int idx) const
    {
        if (idx >= n_local())
        {
            error::out_of_range_error(idx, __FILE__, __LINE__);
        }
        return local_to_global_[idx];
    }
    //=============================================================================
    std::vector<int> IndexMap::local_to_global(const std::vector<int> &idxs) const
    {
        std::vector<int> idxs_mapped(idxs.size());
        for (std::size_t i = 0; i < idxs.size(); i++)
        {
            idxs_mapped[i] = local_to_global(idxs[i]);
        }
        return idxs_mapped;
    }
    //=============================================================================
    int IndexMap::global_to_local(int idx) const
    {
        if (global_to_local_.count(idx) > 0)
        {
            return global_to_local_.at(idx);
        }
        else
        {
            return -1;
        }
    }
    //=============================================================================
    std::vector<int> IndexMap::global_to_local(const std::vector<int> &idxs) const
    {
        std::vector<int> idxs_mapped(idxs.size());
        for (std::size_t i = 0; i < idxs.size(); i++)
        {
            idxs_mapped[i] = global_to_local(idxs[i]);
        }
        return idxs_mapped;
    }
    //=============================================================================
    IndexMap IndexMap::renumber() const
    {
        int proc_rank = Logger::instance().proc_rank();
        int n_procs = Logger::instance().n_procs();

        // If serial, return early
        if (n_procs == 1)
        {
            return IndexMap(n_owned_);
        }

        std::vector<int> owned_idxs_re(n_owned_);
        std::vector<int> ghost_idxs_re(n_ghost_);

        // Renumber the owned indices
        {
            // Compute the displacement for each process
            std::vector<int> send_buffer(n_procs, n_owned_);
            std::vector<int> recv_buffer(n_procs);
            MPI_Alltoall(send_buffer.data(), 1, MPI_INT, recv_buffer.data(), 1, MPI_INT, SFEM_COMM_WORLD);
            int disp = 0;
            for (int i = 0; i < proc_rank; i++)
            {
                disp += recv_buffer[i];
            }

            // Perform the renumbering
            for (int i = 0; i < n_owned_; i++)
            {
                owned_idxs_re[i] = i + disp;
            }
        }

        // Renumber the ghost indices
        {
            int send_buffer_size = 0;
            int recv_buffer_size = 0;
            std::vector<int> send_buffer;
            std::vector<int> recv_buffer(n_ghost_);
            std::vector<int> send_counts(n_procs, 0);
            std::vector<int> recv_counts(n_procs, 0);
            std::vector<int> send_displs(n_procs, 0);
            std::vector<int> recv_displs(n_procs, 0);

            // Compute the send buffer counts and displacements
            for (int i = 0; i < n_ghost_; i++)
            {
                send_counts[ghost_owners_[i]] += 1;
            }
            for (int i = 0; i < n_procs; i++)
            {
                if (i > 0)
                {
                    send_displs[i] = send_displs[i - 1] + send_counts[i - 1];
                }
                send_buffer_size += send_counts[i];
            }
            send_buffer.resize(send_buffer_size);

            // Fill the send buffer
            std::vector<int> offsets(n_procs, 0);
            std::vector<int> position_map(n_ghost_);
            for (int i = 0; i < n_ghost_; i++)
            {
                int owner = ghost_owners_[i];
                send_buffer[send_displs[owner] + offsets[owner]] = local_to_global_[i + n_owned_];
                position_map[i] = send_displs[owner] + offsets[owner];
                offsets[owner]++;
            }

            // Compute the receive buffer size, e.g. the total number of ghost nodes this proces will receive
            MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, SFEM_COMM_WORLD);
            for (int i = 0; i < n_procs; i++)
            {
                if (i > 0)
                {
                    recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
                }
                recv_buffer_size += recv_counts[i];
            }
            recv_buffer.resize(recv_buffer_size);

            // Send the ghost indices to their owner process
            MPI_Alltoallv(send_buffer.data(), send_counts.data(), send_displs.data(), MPI_INT,
                          recv_buffer.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                          SFEM_COMM_WORLD);

            // Renumber the received ghost indices
            for (int i = 0; i < recv_buffer_size; i++)
            {
                int old_global_idx = recv_buffer[i];
                int local_idx = global_to_local_.at(old_global_idx);
                int new_global_idx = owned_idxs_re[local_idx];
                recv_buffer[i] = new_global_idx;
            }

            // Send the renubmered ghost indices back to the process that sent them
            MPI_Alltoallv(recv_buffer.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                          send_buffer.data(), send_counts.data(), send_displs.data(), MPI_INT,
                          SFEM_COMM_WORLD);

            // Correctly place the renumberd ghost indices
            for (int i = 0; i < n_ghost_; i++)
            {
                ghost_idxs_re[i] = send_buffer[position_map[i]];
            }
        }
        return IndexMap(owned_idxs_re, ghost_idxs_re, ghost_owners_);
    }
}