#include "cell.h"
#include "../common/error.h"

namespace sfem::mesh
{
    //=============================================================================
    Cell::Cell(int idx, CellType type, int order, int region_tag)
        : idx_(idx), type_(type), order_(order), region_tag_(region_tag)
    {
        dim_ = cell_dim(type);
        n_nodes_ = cell_nodes(type, order);
        n_faces_ = cell_n_faces(type);
        if (idx < 0 || order < 0 || dim_ < 0 || n_nodes_ < 0 || n_faces_ < 0)
        {
            error::invalid_cell_error(idx, static_cast<int>(type), order, __FILE__, __LINE__);
        }
    }
    //=============================================================================
    int Cell::idx() const
    {
        return idx_;
    }
    //=============================================================================
    CellType Cell::type() const
    {
        return type_;
    }
    //=============================================================================
    int Cell::order() const
    {
        return order_;
    }
    //=============================================================================
    int Cell::region_tag() const
    {
        return region_tag_;
    }
    //=============================================================================
    int Cell::dim() const
    {
        return dim_;
    }
    //=============================================================================
    int Cell::n_nodes() const
    {
        return n_nodes_;
    }
    //=============================================================================
    int Cell::n_faces() const
    {
        return n_faces_;
    }
    //=============================================================================
    geo::Vec3 Cell::face_normal(int f_idx, const std::vector<Scalar> &xpts) const
    {
        // Face nodes
        auto fn = cell_n_nodes_face(type_, f_idx);
        if (fn[0] < 0)
        {
            error::invalid_face_error("", f_idx, __FILE__, __LINE__);
        }

        switch (dim_)
        {
        case 1:
            // Outward normal of line in 2D
            if (f_idx == -1)
            {
                return geo::Vec3(xpts[fn[0] * 3 + 0], xpts[fn[0] * 3 + 1], xpts[fn[0] * 3 + 2],
                                 xpts[fn[1] * 3 + 0], xpts[fn[1] * 3 + 1], xpts[fn[1] * 3 + 2])
                    .unit_norm();
            }
            else
            {
                return geo::Vec3(xpts[fn[0] * 3 + 0], xpts[fn[0] * 3 + 1], xpts[fn[0] * 3 + 2],
                                 xpts[fn[1] * 3 + 0], xpts[fn[1] * 3 + 1], xpts[fn[1] * 3 + 2])
                    .unit_tan();
            }

        case 2:
            // Outward normal of triangle/quad in 3D
            if (f_idx == -1)
            {
                geo::Vec3 v1(xpts[fn[0] * 3 + 0], xpts[fn[0] * 3 + 1], xpts[fn[0] * 3 + 2],
                             xpts[fn[1] * 3 + 0], xpts[fn[1] * 3 + 1], xpts[fn[1] * 3 + 2]);

                geo::Vec3 v2(xpts[fn[2] * 3 + 0], xpts[fn[2] * 3 + 1], xpts[fn[2] * 3 + 2],
                             xpts[fn[3] * 3 + 0], xpts[fn[3] * 3 + 1], xpts[fn[3] * 3 + 2]);

                return v1.cross_prod(v2).normalize();
            }
            else
            {
                return geo::Vec3(xpts[fn[0] * 3 + 0], xpts[fn[0] * 3 + 1], xpts[fn[0] * 3 + 2],
                                 xpts[fn[1] * 3 + 0], xpts[fn[1] * 3 + 1], xpts[fn[1] * 3 + 2])
                    .unit_norm();
            }

        case 3:
            /// @todo Fix
            return geo::Vec3(0, 0, 0);

        default:
            /// @todo Error
            return geo::Vec3(0, 0, 0);
        }
    }
}