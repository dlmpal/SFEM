#pragma once

#include "../common/index_map.h"
#include "../mesh/mesh.h"

using namespace sfem::common;

namespace sfem::la
{
    std::pair<std::vector<int>, std::vector<int>>
    sparsity_pattern(const mesh::Mesh &mesh, int n_vars);
}