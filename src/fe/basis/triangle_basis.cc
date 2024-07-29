#include "basis.h"
#include "quadrature.h"

namespace sfem::fe::basis
{
    //=============================================================================
    int LinearTriangleBasis::dim() const
    {
        return 2;
    }
    //=============================================================================
    int LinearTriangleBasis::n_nodes() const
    {
        return 3;
    }
    //=============================================================================
    int LinearTriangleBasis::n_qpts() const
    {
        return 3;
    }
    //=============================================================================
    Scalar LinearTriangleBasis::qwt(int npt) const
    {
        return triangle_qwt_3[npt];
    }
    //=============================================================================
    void LinearTriangleBasis::qpt(int npt, Scalar pt[]) const
    {
        pt[0] = triangle_qpt_3[npt * 2 + 0];
        pt[1] = triangle_qpt_3[npt * 2 + 1];
    }
    //=============================================================================
    void LinearTriangleBasis::eval_shape(const Scalar pt[], Scalar N[]) const
    {
        N[0] = 1 - pt[0] - pt[1];
        N[1] = pt[0];
        N[2] = pt[1];
    }
    //=============================================================================
    void LinearTriangleBasis::eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const
    {
        dNdxi[0 * 3 + 0] = -1.0;
        dNdxi[0 * 3 + 1] = -1.0;

        dNdxi[1 * 3 + 0] = 1.0;
        dNdxi[1 * 3 + 1] = 0.0;

        dNdxi[2 * 3 + 0] = 0.0;
        dNdxi[2 * 3 + 1] = 1.0;
    }
    //=============================================================================
    int QuadraticTriangleBasis::dim() const
    {
        return 2;
    }
    //=============================================================================
    int QuadraticTriangleBasis::n_nodes() const
    {
        return 6;
    }
    //=============================================================================
    int QuadraticTriangleBasis::n_qpts() const
    {
        return 4;
    }
    //=============================================================================
    Scalar QuadraticTriangleBasis::qwt(int npt) const
    {
        return triangle_qwt_4[npt];
    }
    //=============================================================================
    void QuadraticTriangleBasis::qpt(int npt, Scalar pt[]) const
    {
        pt[0] = triangle_qpt_4[npt * 2 + 0];
        pt[1] = triangle_qpt_4[npt * 2 + 1];
    }
    //=============================================================================
    void QuadraticTriangleBasis::eval_shape(const Scalar pt[], Scalar N[]) const
    {
        // Corner nodes
        N[0] = (1 - pt[0] - pt[1]) * (1 - 2 * pt[0] - 2 * pt[1]);
        N[1] = pt[0] * (2 * pt[0] - 1);
        N[2] = pt[1] * (2 * pt[1] - 1);

        // Mid-side nodes
        N[3] = 4 * pt[0] * (1 - pt[0] - pt[1]);
        N[4] = 4 * pt[0] * pt[1];
        N[5] = 4 * pt[1] * (1 - pt[0] - pt[1]);
    }
    //=============================================================================
    void QuadraticTriangleBasis::eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const
    {
        // Corner nodes
        dNdxi[0 * 3 + 0] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;
        dNdxi[0 * 3 + 1] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;

        dNdxi[1 * 3 + 0] = 4.0 * pt[0] - 1.0;
        dNdxi[1 * 3 + 1] = 0.0;

        dNdxi[2 * 3 + 0] = 0.0;
        dNdxi[2 * 3 + 1] = 4.0 * pt[1] - 1.0;

        // Mid-side nodes
        dNdxi[3 * 3 + 0] = 4.0 - 8.0 * pt[0] - 4.0 * pt[1];
        dNdxi[3 * 3 + 1] = -4.0 * pt[0];

        dNdxi[4 * 3 + 0] = 4.0 * pt[1];
        dNdxi[4 * 3 + 1] = 4.0 * pt[0];

        dNdxi[5 * 3 + 0] = -4.0 * pt[1];
        dNdxi[5 * 3 + 1] = 4.0 - 4.0 * pt[0] - 8.0 * pt[1];
    }
    //=============================================================================
    int CubicTriangleBasis::dim() const
    {
        return 2;
    }
    //=============================================================================
    int CubicTriangleBasis::n_nodes() const
    {
        return 10;
    }
    //=============================================================================
    int CubicTriangleBasis::n_qpts() const
    {
        return 6;
    }
    //=============================================================================
    Scalar CubicTriangleBasis::qwt(int npt) const
    {
        return triangle_qwt_6[npt];
    }
    //=============================================================================
    void CubicTriangleBasis::qpt(int npt, Scalar pt[]) const
    {
        pt[0] = triangle_qpt_6[npt * 2];
        pt[1] = triangle_qpt_6[npt * 2 + 1];
    }
    //=============================================================================
    void CubicTriangleBasis::eval_shape(const Scalar pt[], Scalar N[]) const
    {
        Scalar L1 = 1 - pt[0] - pt[1];
        Scalar L2 = pt[0];
        Scalar L3 = pt[1];

        // Corner nodes
        N[0] = 0.5 * (3 * L1 - 1) * (3 * L1 - 2) * L1;
        N[1] = 0.5 * (3 * L2 - 1) * (3 * L2 - 2) * L2;
        N[2] = 0.5 * (3 * L3 - 1) * (3 * L3 - 2) * L3;

        // Mid-side nodes
        N[3] = 9 / 2 * L1 * L2 * (3 * L1 - 1);
        N[4] = 9 / 2 * L1 * L2 * (3 * L2 - 1);

        N[5] = 9 / 2 * L2 * L3 * (3 * L2 - 1);
        N[6] = 9 / 2 * L2 * L3 * (3 * L3 - 1);

        N[7] = 9 / 2 * L3 * L1 * (3 * L3 - 1);
        N[8] = 9 / 2 * L3 * L1 * (3 * L1 - 1);

        // Internal node
        N[9] = 27 * L1 * L2 * L3;
    }
    //=============================================================================
    void CubicTriangleBasis::eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const
    {
        // Corner nodes
        dNdxi[0] = -13.5 * pt[1] * pt[1] - 27.0 * pt[1] * pt[0] + 18.0 * pt[1] - 13.5 * pt[0] * pt[0] + 18.0 * pt[0] - 5.5;
        dNdxi[1] = -13.5 * pt[1] * pt[1] - 27.0 * pt[1] * pt[0] + 18.0 * pt[1] - 13.5 * pt[0] * pt[0] + 18.0 * pt[0] - 5.5;

        dNdxi[3] = 13.5 * pt[0] * pt[0] - 9.0 * pt[0] + 1.0;
        dNdxi[4] = 0;

        dNdxi[6] = 0;
        dNdxi[7] = 13.5 * pt[1] * pt[1] - 9.0 * pt[1] + 1.0;

        // Mid-side nodes
        dNdxi[9] = 13.5 * pt[1] * pt[1] + 54.0 * pt[1] * pt[0] - 22.5 * pt[1] + 40.5 * pt[0] * pt[0] - 45.0 * pt[0] + 9.0;
        dNdxi[10] = pt[0] * (27.0 * pt[1] + 27.0 * pt[0] - 22.5);

        dNdxi[12] = -27.0 * pt[1] * pt[0] + 4.5 * pt[1] - 40.5 * pt[0] * pt[0] + 36.0 * pt[0] - 4.5;
        dNdxi[13] = pt[0] * (4.5 - 13.5 * pt[0]);

        dNdxi[15] = pt[1] * (27.0 * pt[0] - 4.5);
        dNdxi[16] = pt[0] * (13.5 * pt[0] - 4.5);

        dNdxi[18] = pt[1] * (13.5 * pt[1] - 4.5);
        dNdxi[19] = pt[0] * (27.0 * pt[1] - 4.5);

        dNdxi[21] = pt[1] * (4.5 - 13.5 * pt[1]);
        dNdxi[22] = -40.5 * pt[1] * pt[1] - 27.0 * pt[1] * pt[0] + 36.0 * pt[1] + 4.5 * pt[0] - 4.5;

        dNdxi[24] = pt[1] * (27.0 * pt[1] + 27.0 * pt[0] - 22.5);
        dNdxi[25] = 40.5 * pt[1] * pt[1] + 54.0 * pt[1] * pt[0] - 45.0 * pt[1] + 13.5 * pt[0] * pt[0] - 22.5 * pt[0] + 9.0;

        // Internal node
        dNdxi[27] = 27 * pt[1] * (-pt[1] - 2 * pt[0] + 1);
        dNdxi[28] = 27 * pt[0] * (-2 * pt[1] - pt[0] + 1);
    }
}