#include "basis.h"
#include "quadrature.h"

namespace sfem::fe::basis
{
    //=============================================================================
    int LinearHexahedralBasis::dim() const
    {
        return 3;
    }
    //=============================================================================
    int LinearHexahedralBasis::n_nodes() const
    {
        return 8;
    }
    //=============================================================================
    int LinearHexahedralBasis::n_qpts() const
    {
        return 8;
    }
    //=============================================================================
    Scalar LinearHexahedralBasis::qwt(int npt) const
    {
        return gauss_qwt_2[npt % 2] * gauss_qwt_2[(npt % 4) / 2] * gauss_qwt_2[npt / 4];
    }
    //=============================================================================
    void LinearHexahedralBasis::qpt(int npt, Scalar pt[]) const
    {
        pt[0] = gauss_qpt_2[npt % 2];
        pt[1] = gauss_qpt_2[(npt % 4) / 2];
        pt[2] = gauss_qpt_2[npt / 2];
    }
    //=============================================================================
    void LinearHexahedralBasis::eval_shape(const Scalar pt[], Scalar N[]) const
    {
        N[0] = (0.125 - 0.125 * pt[0]) * (1 - pt[1]) * (1 - pt[2]);
        N[1] = (1 - pt[1]) * (1 - pt[2]) * (0.125 * pt[0] + 0.125);
        N[2] = (1 - pt[2]) * (pt[1] + 1) * (0.125 * pt[0] + 0.125);
        N[3] = (0.125 - 0.125 * pt[0]) * (1 - pt[2]) * (pt[1] + 1);
        N[4] = (0.125 - 0.125 * pt[0]) * (1 - pt[1]) * (pt[2] + 1);
        N[5] = (1 - pt[1]) * (0.125 * pt[0] + 0.125) * (pt[2] + 1);
        N[6] = (pt[1] + 1) * (0.125 * pt[0] + 0.125) * (pt[2] + 1);
        N[7] = (0.125 - 0.125 * pt[0]) * (pt[1] + 1) * (pt[2] + 1);
    }
    //=============================================================================
    void LinearHexahedralBasis::eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const
    {
        dNdxi[0] = -0.125 * (pt[1] - 1) * (pt[2] - 1);
        dNdxi[1] = -0.125 * (pt[0] - 1) * (pt[2] - 1);
        dNdxi[2] = -0.125 * (pt[1] - 1) * (pt[0] - 1);

        dNdxi[3] = 0.125 * (pt[1] - 1) * (pt[2] - 1);
        dNdxi[4] = 0.125 * (pt[0] + 1) * (pt[2] - 1);
        dNdxi[5] = 0.125 * (pt[1] - 1) * (pt[0] + 1);

        dNdxi[6] = -0.125 * (pt[1] + 1) * (pt[2] - 1);
        dNdxi[7] = -0.125 * (pt[0] + 1) * (pt[2] - 1);
        dNdxi[8] = -0.125 * (pt[1] + 1) * (pt[0] + 1);

        dNdxi[9] = 0.125 * (pt[1] + 1) * (pt[2] - 1);
        dNdxi[10] = 0.125 * (pt[0] - 1) * (pt[2] - 1);
        dNdxi[11] = 0.125 * (pt[1] + 1) * (pt[0] - 1);

        dNdxi[12] = 0.125 * (pt[1] - 1) * (pt[2] + 1);
        dNdxi[13] = 0.125 * (pt[0] - 1) * (pt[2] + 1);
        dNdxi[14] = 0.125 * (pt[1] - 1) * (pt[0] - 1);

        dNdxi[15] = -0.125 * (pt[1] - 1) * (pt[2] + 1);
        dNdxi[16] = -0.125 * (pt[0] + 1) * (pt[2] + 1);
        dNdxi[17] = -0.125 * (pt[1] - 1) * (pt[0] + 1);

        dNdxi[18] = 0.125 * (pt[1] + 1) * (pt[2] + 1);
        dNdxi[19] = 0.125 * (pt[0] + 1) * (pt[2] + 1);
        dNdxi[20] = 0.125 * (pt[1] + 1) * (pt[0] + 1);

        dNdxi[21] = -0.125 * (pt[1] + 1) * (pt[2] + 1);
        dNdxi[22] = -0.125 * (pt[0] - 1) * (pt[2] + 1);
        dNdxi[23] = -0.125 * (pt[1] + 1) * (pt[0] - 1);
    }
}