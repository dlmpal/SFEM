#include "mpi_utils.h"
#include "logger.h"
#include "error.h"
#include <mpi.h>

namespace sfem::mpi
{
    //=============================================================================
    std::vector<int> send_data_to_owners(const std::vector<int> owners, const std::vector<int> data)
    {
        if (owners.size() != data.size())
        {
            error::invalid_size_error(owners.size(), data.size(), __FILE__, __LINE__);
        }

        int n_procs = Logger::instance().n_procs();
        // If serial, return early
        if (n_procs == 1)
        {
            std::vector<int>(0);
        }

        // Compute the send counts, e.g.
        // home many entries this process sends to each other
        std::vector<int> send_counts(n_procs, 0);
        for (std::size_t i = 0; i < owners.size(); i++)
        {
            send_counts[owners[i]]++;
        }

        // Compute the send displacements
        std::vector<int> send_displs(n_procs, 0);
        for (int i = 1; i < n_procs; i++)
        {
            send_displs[i] = send_displs[i - 1] + send_counts[i - 1];
        }

        // Allocate and fill the send buffer
        std::vector<int> send_buffer(data.size());
        std::vector<int> offsets(n_procs, 0);
        for (std::size_t i = 0; i < data.size(); i++)
        {
            send_buffer[send_displs[owners[i]] + offsets[owners[i]]] = data[i];
            offsets[owners[i]]++;
        }

        // Compute the receive counts, e.g.
        // how many entries this process receives from each other
        std::vector<int> recv_counts(n_procs, 0);
        MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                     recv_counts.data(), 1, MPI_INT, SFEM_COMM_WORLD);

        // Compute the receive displacements
        std::vector<int> recv_displs(n_procs, 0);
        for (int i = 1; i < n_procs; i++)
        {
            recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
        }

        // Allocate the receive buffer
        int recv_buffer_size = 0;
        for (int i = 0; i < n_procs; i++)
        {
            recv_buffer_size += recv_counts[i];
        }
        std::vector<int> recv_buffer(recv_buffer_size);

        // Send and receive data to/from all processes
        MPI_Alltoallv(send_buffer.data(), send_counts.data(), send_displs.data(), MPI_INT,
                      recv_buffer.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                      SFEM_COMM_WORLD);

        return recv_buffer;
    }
}