#pragma once

#include "../mesh/field.h"

namespace sfem::io
{
    /// @brief Read field values from file
    void read_field_values(const std::string &path, mesh::Field &field);

    /// @brief Write field values to file
    void write_field_values(const std::string &path, const mesh::Field &field, bool assemble_global = true);
}