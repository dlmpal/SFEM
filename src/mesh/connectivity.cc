#include "connectivity.h"

namespace sfem::mesh
{
    //=============================================================================
    Connectivity invert_conn(const Connectivity &conn)
    {
        Connectivity inv;
        inv.n1 = conn.n2;
        inv.n2 = conn.n1;

        inv.cnt.resize(inv.n1);
        std::fill(inv.cnt.begin(), inv.cnt.end(), 0);
        for (int i = 0; i < conn.n1; i++)
        {
            for (int j = 0; j < conn.cnt[i]; j++)
            {
                int idx = conn.idx.at(conn.ptr.at(i) + j);
                inv.cnt[idx]++;
            }
        }

        inv.ptr.resize(inv.n1);
        int ptr = 0;
        for (int i = 0; i < inv.n1; i++)
        {
            inv.ptr[i] = ptr;
            ptr += inv.cnt[i];
        }

        std::fill(inv.cnt.begin(), inv.cnt.end(), 0);
        inv.idx.resize(conn.idx.size());
        for (int i = 0; i < conn.n1; i++)
        {
            for (int j = 0; j < conn.cnt[i]; j++)
            {
                int idx = conn.idx.at(conn.ptr.at(i) + j);
                inv.idx[inv.ptr.at(idx) + inv.cnt.at(idx)] = i;
                inv.cnt[idx]++;
            }
        }

        return inv;
    }
    //=============================================================================
    Connectivity compute_node_to_node_conn(const Connectivity &cell_node_conn)
    {
        Connectivity node_cell_conn = invert_conn(cell_node_conn);

        Connectivity node_node_conn;
        node_node_conn.n1 = node_cell_conn.n1;
        node_node_conn.n2 = node_cell_conn.n1;
        node_node_conn.cnt.resize(node_node_conn.n1);
        node_node_conn.ptr.resize(node_node_conn.n1);

        // First loop to get count
        // e.g. how many nodes are surrounding each node
        {
            // Loop over the nodes
            int size = 0;
            std::vector<int> included(node_node_conn.n1, -1);
            for (int i = 0; i < node_cell_conn.n1; i++)
            {
                // Loop over the cells surrounding the node
                for (int j = 0; j < node_cell_conn.cnt[i]; j++)
                {
                    int cell_idx = node_cell_conn.idx.at(node_cell_conn.ptr.at(i) + j);

                    // Loop over the nodes of the cell
                    for (int k = 0; k < cell_node_conn.cnt.at(cell_idx); k++)
                    {
                        int node_idx = cell_node_conn.idx.at(cell_node_conn.ptr.at(cell_idx) + k);

                        if (included.at(node_idx) != i)
                        {
                            size++;
                            node_node_conn.cnt[i]++;
                            included[node_idx] = i;
                        }
                    }
                }
            }
            node_node_conn.idx.resize(size);
        }

        int ptr = 0;
        for (int i = 0; i < node_node_conn.n1; i++)
        {
            node_node_conn.ptr[i] = ptr;
            ptr += node_node_conn.cnt[i];
        }

        // Second loop to assign nodes
        {
            // Loop over the nodes
            std::fill(node_node_conn.cnt.begin(), node_node_conn.cnt.end(), 0);
            std::vector<int> included(node_cell_conn.n1, -1);
            for (int i = 0; i < node_cell_conn.n1; i++)
            {
                // Loop over the cells surrounding the node
                for (int j = 0; j < node_cell_conn.cnt[i]; j++)
                {
                    int cell_idx = node_cell_conn.idx.at(node_cell_conn.ptr.at(i) + j);

                    // Loop over the nodes of the cell
                    for (int k = 0; k < cell_node_conn.cnt.at(cell_idx); k++)
                    {
                        int node_idx = cell_node_conn.idx.at(cell_node_conn.ptr.at(cell_idx) + k);

                        if (included.at(node_idx) != i)
                        {
                            node_node_conn.idx[node_node_conn.ptr.at(i) + node_node_conn.cnt.at(i)] = node_idx;
                            node_node_conn.cnt[i]++;
                            included[node_idx] = i;
                        }
                    }
                }
            }
        }

        return node_node_conn;
    }
}