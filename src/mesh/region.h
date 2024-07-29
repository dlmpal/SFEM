#pragma once

#include "../common/config.h"
#include <string>

namespace sfem::mesh
{
  class Region
  {
  public:
    /// @brief Create a region
    Region(const std::string &name, int dim, int tag);

    std::string name() const;
    int dim() const;
    int tag() const;

  private:
    /// @brief Name
    std::string name_;

    /// @brief Physical dimension
    int dim_;

    /// @brief Numerical tag
    int tag_;
  };
}
