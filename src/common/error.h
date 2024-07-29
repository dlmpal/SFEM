#pragma once

#include "logger.h"

/// @brief Commonly used errors
namespace sfem::error
{
    void invalid_filename_error(const std::string &filename, const std::string &file, int line);
    void invalid_size_error(int correct_size, int wrong_size, const std::string &file, int line);
    void out_of_range_error(int idx, const std::string &file, int line);
    void invalid_cell_error(int cell_id, int cell_type, int cell_order, const std::string &file, int line);
    void negative_jacobian_error(int cell_id, const std::string &file, int line);
    void invalid_face_error(const std::string cell_type, int f_idx, const std::string &file, int line);
    void unsupported_gmsh_type_error(int gmsh_id, int gmsh_type, const std::string &file, int line);
}