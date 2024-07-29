#include "error.h"

namespace sfem::error
{
    //=============================================================================
    void invalid_filename_error(const std::string &filename, const std::string &file, int line)
    {
        std::string message = "Error opening file: " + filename + "\n";
        Logger::instance().error(message, file, line);
    }
    //=============================================================================
    void invalid_size_error(int correct_size, int wrong_size, const std::string &file, int line)
    {
        std::string message = "Expected size: " + std::to_string(correct_size) + "\n";
        message += "Got size: " + std::to_string(wrong_size) + "\n";
        Logger::instance().error(message, file, line);
    }
    //=============================================================================
    void out_of_range_error(int idx, const std::string &file, int line)
    {
        std::string message = "Index: " + std::to_string(idx) + " is out of range\n";
        Logger::instance().error(message, file, line);
    }
    //=============================================================================
    void invalid_cell_error(int cell_id, int cell_type, int cell_order, const std::string &file, int line)
    {
        std::string message = "Cell " + std::to_string(cell_id) + " has invalid type (" + std::to_string(cell_type) + ") or order (" + std::to_string(cell_order) + ")\n";
        Logger::instance().error(message, file, line);
    }
    //=============================================================================
    void negative_jacobian_error(int cell_id, const std::string &file, int line)
    {
        std::string message = "Negative jacobian at cell: " + std::to_string(cell_id) + "\n";
        Logger::instance().error(message, file, line);
    }
    //=============================================================================
    void invalid_face_error(std::string cell_type, int f_idx, const std::string &file, int line)
    {
        std::string message = "Invalid face index: " + std::to_string(f_idx) + "for " + cell_type + "\n";
        Logger::instance().error(message, file, line);
    }
    //=============================================================================
    void unsupported_gmsh_type_error(int gmsh_id, int gmsh_type, const std::string &file, int line)
    {
        std::string message = "Unsupported Gmsh cell type: " + std::to_string(gmsh_type) + "\n";
        Logger::instance().error(message, file, line);
    }
}