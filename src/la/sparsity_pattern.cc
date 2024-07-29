#include "sparsity_pattern.h"
#include "../common/mpi_utils.h"
#include "../common/timer.h"
#include <numeric>

namespace sfem::la
{
    //=============================================================================
    std::pair<std::vector<int>, std::vector<int>> sparsity_pattern(const mesh::Mesh &mesh, int n_vars)
    {
        Timer timer("Sparsity pattern calculation");

        // Node index map (renumbered)
        auto im = mesh.node_im().renumber();

        // Node-to-node connectivity
        auto conn = mesh::compute_node_to_node_conn(mesh.cell_node_conn());

        // Number of non-zeros for locally owned indices
        std::vector<int> diag_nnz(im.n_owned(), 0);
        std::vector<int> off_diag_nnz(im.n_owned(), 0);

        // Number of non-zeros for ghost indices
        // These have to be sent to the ghost indices' owners
        std::vector<int> ghost_diag_nnz(im.n_ghost(), 0);
        std::vector<int> ghost_off_diag_nnz(im.n_ghost(), 0);

        // Get the ghost indices and their owners
        auto ghost_idxs = im.get_ghost_idxs();
        auto ghost_owners = im.get_ghost_owners();

        // Loop over all indices
        for (int i = 0; i < conn.n1; i++)
        {
            // Loop over every neighbouring index
            for (int k = 0; k < conn.cnt[i]; k++)
            {
                int j = conn.idx[conn.ptr[i] + k];

                // Index "i" is locally owned
                if (i < im.n_owned())
                {
                    // Index "j" is locally owned
                    if (j < im.n_owned())
                    {
                        diag_nnz[i]++;
                    }
                    else
                    {
                        off_diag_nnz[i]++;
                    }
                }
                // Index "i" is owned by another process
                else
                {
                    int owner_i = ghost_owners[i - im.n_owned()];

                    // Get the owner of index "j"
                    // owner_j = -1 if "j" is locally owned
                    int owner_j;
                    if (j < im.n_owned())
                    {
                        owner_j = -1;
                    }
                    else
                    {
                        owner_j = ghost_owners[j - im.n_owned()];
                    }

                    // Indices "i" and "j" are owned by the same process
                    if (owner_i == owner_j)
                    {
                        ghost_diag_nnz[i - im.n_owned()]++;
                    }
                    else
                    {
                        ghost_off_diag_nnz[i - im.n_owned()]++;
                    }
                }
            }
        }

        // Locally owned indices, that are ghosts for other processes
        auto a1 = mpi::send_data_to_owners(ghost_owners, ghost_idxs);

        // diag and off-diag non-zeros computed for locally owned indices
        // by other processes
        auto a2 = mpi::send_data_to_owners(ghost_owners, ghost_diag_nnz);
        auto a3 = mpi::send_data_to_owners(ghost_owners, ghost_off_diag_nnz);

        // Add the non-zeros computed non-locally to the local ones
        for (std::size_t i = 0; i < a1.size(); i++)
        {
            diag_nnz[im.global_to_local(a1[i])] += a2[i];
            off_diag_nnz[im.global_to_local(a1[i])] += a3[i];
        }

        // Resize the NNZ vectors for the given number of variables per index
        diag_nnz.resize(im.n_owned() * n_vars);
        off_diag_nnz.resize(im.n_owned() * n_vars);
        for (int i = im.n_owned() - 1; i >= 0; i--)
        {
            for (int j = 0; j < n_vars; j++)
            {
                diag_nnz[i * n_vars + j] = diag_nnz[i] * n_vars;
                off_diag_nnz[i * n_vars + j] = off_diag_nnz[i] * n_vars;
            }
        }

        // Print the number of non-zeros
        int n_diag_nnz = 0;
        int n_off_diag_nnz = 0;
        n_diag_nnz = std::accumulate(diag_nnz.begin(), diag_nnz.end(), n_diag_nnz);
        n_off_diag_nnz = std::accumulate(off_diag_nnz.begin(), off_diag_nnz.end(), n_off_diag_nnz);
        std::string message = "Number of non-zeros, Diag: " + std::to_string(n_diag_nnz) + " Off-Diag: " + std::to_string(n_off_diag_nnz) + "\n";
        Logger::instance().log_message(message, Logger::Level::all);

        return std::make_pair(diag_nnz, off_diag_nnz);
    }
}