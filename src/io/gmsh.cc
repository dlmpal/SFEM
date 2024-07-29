#include "gmsh.h"
#include <algorithm>

namespace sfem::io
{
    //=============================================================================
    mesh::Mesh read_gmsh(const std::string &path)
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            error::invalid_filename_error(path, __FILE__, __LINE__);
        }

        // Keep track of current file line
        int line = 0;

        // Skip first 4 lines
        std::string buffer;
        for (int i = 0; i < 6; i++)
        {
            file >> buffer;
            line++;
        }

        // Read regions ($PhysicalNames)
        int n_regions;
        std::vector<mesh::Region> regions;
        file >> n_regions;
        line++;
        for (int i = 0; i < n_regions; i++)
        {
            std::string name;
            int dim, tag;
            file >> dim;
            file >> tag;
            file >> name;
            name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
            regions.push_back(mesh::Region(name, dim, tag));
            line += 3;
        }

        // Skip two lines
        for (int i = 0; i < 2; i++)
        {
            file >> buffer;
            line++;
        }

        // Read nodal positions ($Nodes)
        int n_nodes;
        file >> n_nodes;
        line++;
        std::vector<Scalar> xpts(n_nodes * 3);
        for (int i = 0; i < n_nodes; i++)
        {
            file >> buffer;
            for (int j = 0; j < 3; j++)
            {
                file >> xpts[i * 3 + j];
            }
            line += 4;
        }

        // Skip two lines
        for (int i = 0; i < 2; i++)
        {
            file >> buffer;
            line++;
        }

        // Read elements ($Elements), but skip the connectivity
        int n_cells, size_conn = 0;
        file >> n_cells;
        line++;
        std::vector<mesh::Cell> cells;
        cells.reserve(n_cells);
        for (int i = 0; i < n_cells; i++)
        {
            int id, gmsh_type, dim, physical_tag, elementary_tag;

            file >> id;
            file >> gmsh_type;
            file >> dim;
            file >> physical_tag;
            file >> elementary_tag;
            auto [type, order] = gmsh::gmsh_type_to_native(id, gmsh_type);
            mesh::Cell cell(id - 1, type, order, physical_tag);
            for (int j = 0; j < cell.n_nodes(); j++)
            {
                file >> buffer;
            }
            size_conn += cell.n_nodes();
            cells.push_back(cell);
        }

        // Close the file and re-open it
        // at the $Elements section
        file.close();
        file.open(path);
        for (int i = 0; i < line; i++)
        {
            file >> buffer;
        }

        // Read the cell-to-node connectivity
        mesh::Connectivity conn;
        conn.n1 = n_cells;
        conn.n2 = n_nodes;
        conn.idx.resize(size_conn);
        conn.ptr.resize(n_cells);
        conn.cnt.resize(n_cells);
        int ptr = 0;
        for (int i = 0; i < n_cells; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                file >> buffer;
            }

            conn.ptr[i] = ptr;
            conn.cnt[i] = cells[i].n_nodes();

            for (int j = 0; j < cells[i].n_nodes(); j++)
            {
                file >> conn.idx[ptr + j];
            }

            ptr += cells[i].n_nodes();
        }

        // Gmsh numbering starts at 1
        for (int i = 0; i < size_conn; i++)
        {
            conn.idx[i] = conn.idx[i] - 1;
        }

        return mesh::Mesh(cells, conn, xpts, regions, common::IndexMap(n_cells), common::IndexMap(n_nodes));
    }
}