#include "vtk.h"

namespace sfem::io
{
    //=============================================================================
    void write_vtk(const std::string &path, const mesh::Mesh &mesh, const std::vector<mesh::Field> &fields)
    {
        std::ofstream file(path);
        if (!file.is_open())
        {
            error::invalid_filename_error(path, __FILE__, __LINE__);
        }

        // File header
        file << "# vtk DataFile Version 2.0\n";
        file << "SFEM\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n";

        // Points
        const auto &xpts = mesh.xpts();
        file << "POINTS " << xpts.size() / 3 << " float\n";
        for (std::size_t i = 0; i < xpts.size() / 3; i++)
        {
            file << xpts[i * 3 + 0] << " " << xpts[i * 3 + 1] << " " << xpts[i * 3 + 2] << "\n";
        }

        // Cells
        const auto &cells = mesh.cells();
        const auto &conn = mesh.cell_node_conn();
        int vtk_size = conn.idx.size() + cells.size();
        file << "CELLS " << cells.size() << " " << vtk_size << "\n";
        for (const auto &cell : cells)
        {
            file << cell.n_nodes() << " ";
            auto cell_nodes = mesh.get_cell_nodes(cell);
            vtk::cell_node_ordering_to_vtk(cell, cell_nodes.data());
            for (int j = 0; j < cell.n_nodes(); j++)
            {
                file << cell_nodes[j] << " ";
            }
            file << "\n";
        }

        // Cell types
        file << "CELL_TYPES " << cells.size() << "\n";
        for (const auto &cell : cells)
        {
            file << vtk::cell_type_to_vtk(cell) << "\n";
        }

        // Field values
        if (fields.size() > 0)
        {
            file << "POINT_DATA " << mesh.n_nodes_local() << "\n";
            for (const auto &field : fields)
            {
                int n_vars = field.n_vars();
                auto comp_names = field.comp_names();
                const auto &values = field.values();
                for (int j = 0; j < n_vars; j++)
                {
                    file << "SCALARS " << " " << comp_names[j] << " float\n";
                    file << "LOOKUP_TABLE default\n";
                    for (int i = 0; i < mesh.n_nodes_local(); i++)
                    {
                        file << values[i * n_vars + j] << "\n";
                    }
                }
            }
        }
    }
}