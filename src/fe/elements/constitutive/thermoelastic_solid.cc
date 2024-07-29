#include "thermoelastic_solid.h"
#include <cmath>

namespace sfem::fe::constitutive
{
    //=============================================================================
    ThermoElasticSolidConstitutive::ThermoElasticSolidConstitutive(ThermoMechanicalProperties &prop)
        : ThermoElasticConstitutive(prop)
    {
    }
    //=============================================================================
    la::DenseMatrix ThermoElasticSolidConstitutive::stress_strain_matrix() const
    {
        la::DenseMatrix D(6, 6);
        Scalar c1 = _prop.E / ((1 + _prop.nu) * (1 - 2 * _prop.nu));
        Scalar c2 = (1 - 2 * _prop.nu) / 2;

        D.insert(0, 0, (1 - _prop.nu) * c1);
        D.insert(0, 1, _prop.nu * c1);
        D.insert(0, 2, _prop.nu * c1);

        D.insert(1, 0, _prop.nu * c1);
        D.insert(1, 1, (1 - _prop.nu) * c1);
        D.insert(1, 2, _prop.nu * c1);

        D.insert(2, 0, _prop.nu * c1);
        D.insert(2, 1, _prop.nu * c1);
        D.insert(2, 2, (1 - _prop.nu) * c1);

        D.insert(3, 3, c1 * c2);
        D.insert(4, 4, c1 * c2);
        D.insert(5, 5, c1 * c2);

        return D;
    }
    //=============================================================================
    la::DenseMatrix ThermoElasticSolidConstitutive::strain_displacement_matrix(const std::vector<Scalar> &dNdX) const
    {
        int n_nodes = dNdX.size() / 3;
        int n_cols = n_nodes * 3;
        int n_rows = 6;
        la::DenseMatrix B(n_rows, n_cols);
        for (int i = 0; i < n_nodes; i++)
        {
            // exx
            B.insert(0, i * 3 + 0, dNdX[i * 3 + 0]);

            // eyy
            B.insert(1, i * 3 + 1, dNdX[i * 3 + 1]);

            // ezz
            B.insert(2, i * 3 + 2, dNdX[i * 3 + 2]);

            // exy
            B.insert(3, i * 3 + 0, dNdX[i * 3 + 1]);
            B.insert(3, i * 3 + 1, dNdX[i * 3 + 0]);

            // eyz
            B.insert(4, i * 3 + 1, dNdX[i * 3 + 2]);
            B.insert(4, i * 3 + 2, dNdX[i * 3 + 1]);

            // exz
            B.insert(5, i * 3 + 0, dNdX[i * 3 + 2]);
            B.insert(5, i * 3 + 2, dNdX[i * 3 + 0]);
        }
        return B;
    }
    //=============================================================================
    Scalar ThermoElasticSolidConstitutive::eval_vm_stress(const std::vector<Scalar> &dNdx, const std::vector<Scalar> &u) const
    {
        auto stress = eval_stress(dNdx, u);

        Scalar s_xx = stress.at(0, 0);
        Scalar s_yy = stress.at(1, 0);
        Scalar s_zz = stress.at(2, 0);
        Scalar s_xy = stress.at(3, 0);
        Scalar s_yz = stress.at(4, 0);
        Scalar s_xz = stress.at(5, 0);

        Scalar stress_vm = pow(s_xx - s_yy, 2) + pow(s_yy - s_zz, 2) + pow(s_zz - s_xx, 2);
        stress_vm += 6 * (pow(s_xy, 2) + pow(s_yz, 2) + pow(s_xz, 2));
        stress_vm = sqrt(0.5 * stress_vm);

        return stress_vm;
    }
    //=============================================================================
    la::DenseMatrix ThermoElasticSolidConstitutive::eval_thermal_strain(const std::vector<Scalar> &N, const std::vector<Scalar> &dT) const
    {
        Scalar alpha = _prop.alpha;
        la::DenseMatrix eps(6, 1);
        for (std::size_t i = 0; i < N.size(); i++)
        {
            eps.add(0, 0, alpha * N[i] * dT[i]);
            eps.add(1, 0, alpha * N[i] * dT[i]);
            eps.add(2, 0, alpha * N[i] * dT[i]);
        }
        return eps;
    }
}