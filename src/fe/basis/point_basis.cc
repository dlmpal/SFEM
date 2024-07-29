#include "basis.h"

namespace sfem::fe::basis
{
    //=============================================================================
    int PointBasis::dim() const
    {
        return 0;
    }
    //=============================================================================
    int PointBasis::n_nodes() const
    {
        return 1;
    }
    //=============================================================================
    int PointBasis::n_qpts() const
    {
        return 1;
    }
    //=============================================================================
    Scalar PointBasis::qwt(int npt) const
    {
        return 1.0;
    }
    //=============================================================================
    void PointBasis::qpt(int npt, Scalar pt[]) const
    {
        pt[0] = 1.0;
    }
    //=============================================================================
    void PointBasis::eval_shape(const Scalar pt[], Scalar N[]) const
    {
        N[0] = 1.0;
    }
    //=============================================================================
    void PointBasis::eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const
    {
        dNdxi[0 * 3 + 0] = 1.0;
    }
}
