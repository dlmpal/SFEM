#include "region.h"

namespace sfem::mesh
{
    //=============================================================================
    Region::Region(const std::string &name, int dim, int tag)
        : name_(name), dim_(dim), tag_(tag)
    {
    }
    //=============================================================================
    std::string Region::name() const
    {
        return name_;
    }
    //=============================================================================
    int Region::dim() const
    {
        return dim_;
    }
    //=============================================================================
    int Region::tag() const
    {
        return tag_;
    }
}