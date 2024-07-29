#include "basis.h"
#include "quadrature.h"

namespace sfem::fe::basis
{
    //=============================================================================
    int LinearLineBasis::dim() const
    {
        return 1;
    }
    //=============================================================================
    int LinearLineBasis::n_nodes() const
    {
        return 2;
    }
    //=============================================================================
    int LinearLineBasis::n_qpts() const
    {
        return 2;
    }
    //=============================================================================
    Scalar LinearLineBasis::qwt(int npt) const
    {
        return gauss_qwt_2[npt];
    }
    //=============================================================================
    void LinearLineBasis::qpt(int npt, Scalar pt[]) const
    {
        pt[0] = gauss_qpt_2[npt];
    }
    //=============================================================================
    void LinearLineBasis::eval_shape(const Scalar pt[], Scalar N[]) const
    {
        N[0] = 0.5 * (1.0 - pt[0]);
        N[1] = 0.5 * (1.0 + pt[0]);
    }
    //=============================================================================
    void LinearLineBasis::eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const
    {
        dNdxi[0 * 3 + 0] = -0.5;
        dNdxi[1 * 3 + 0] = 0.5;
    }
    //=============================================================================
    int QuadraticLineBasis::dim() const
    {
        return 1;
    }
    //=============================================================================
    int QuadraticLineBasis::n_nodes() const
    {
        return 3;
    }
    //=============================================================================
    int QuadraticLineBasis::n_qpts() const
    {
        return 3;
    }
    //=============================================================================
    Scalar QuadraticLineBasis::qwt(int npt) const
    {
        return gauss_qwt_3[npt];
    }
    //=============================================================================
    void QuadraticLineBasis::qpt(int npt, Scalar pt[]) const
    {
        pt[0] = gauss_qpt_3[npt];
    }
    //=============================================================================
    void QuadraticLineBasis::eval_shape(const Scalar pt[], Scalar N[]) const
    {
        N[0] = -0.5 * pt[0] * (1.0 - pt[0]);
        N[1] = 0.5 * pt[0] * (1 + pt[0]);
        N[2] = 1.0 - pt[0] * pt[0];
    }
    //=============================================================================
    void QuadraticLineBasis::eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const
    {
        dNdxi[0 * 3 + 0] = -0.5 + pt[0];
        dNdxi[1 * 3 + 0] = 0.5 + pt[0];
        dNdxi[2 * 3 + 0] = -2.0 * pt[0];
    }
    //=============================================================================
    int CubicLineBasis::dim() const
    {
        return 1;
    }
    //=============================================================================
    int CubicLineBasis::n_nodes() const
    {
        return 4;
    }
    //=============================================================================
    int CubicLineBasis::n_qpts() const
    {
        return 4;
    }
    //=============================================================================
    Scalar CubicLineBasis::qwt(int npt) const
    {
        return gauss_qwt_4[npt];
    }
    //=============================================================================
    void CubicLineBasis::qpt(int npt, Scalar pt[]) const
    {
        pt[0] = gauss_qpt_4[npt];
    }
    //=============================================================================
    void CubicLineBasis::eval_shape(const Scalar pt[], Scalar N[]) const
    {
        N[0] = (0.333333333333333 - pt[0]) * (0.5625 * pt[0] - 0.5625) * (pt[0] + 0.333333333333333);
        N[1] = (0.333333333333333 - pt[0]) * (-0.5625 * pt[0] - 0.5625) * (pt[0] + 0.333333333333333);
        N[2] = (0.333333333333333 - pt[0]) * (1.6875 - 1.6875 * pt[0]) * (pt[0] + 1);
        N[3] = (1.6875 - 1.6875 * pt[0]) * (pt[0] + 0.333333333333333) * (pt[0] + 1);
    }
    //=============================================================================
    void CubicLineBasis::eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const
    {
        dNdxi[0] = -1.6875 * pt[0] * pt[0] + 1.125 * pt[0] + 0.0625;
        dNdxi[3] = 1.6875 * pt[0] * pt[0] + 1.125 * pt[0] - 0.0625;
        dNdxi[6] = 5.0625 * pt[0] * pt[0] - 1.125 * pt[0] - 1.6875;
        dNdxi[9] = -5.0625 * pt[0] * pt[0] - 1.125 * pt[0] + 1.6875;
    }
}