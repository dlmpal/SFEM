#pragma once

#include <vector>

namespace sfem::mpi
{
    std::vector<int> send_data_to_owners(const std::vector<int> owners, const std::vector<int> data);
}