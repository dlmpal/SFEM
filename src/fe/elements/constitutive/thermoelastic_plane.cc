#include "thermoelastic_plane.h"

namespace sfem::fe::constitutive
{
    //=============================================================================
    ThermoElasticPlaneConstitutive::ThermoElasticPlaneConstitutive(
        ThermoMechanicalProperties &prop, Scalar thick, ThermoElasticPlaneConstitutive::Type type)
        : ThermoElasticConstitutive(prop), _thick(thick), _type(type)
    {
    }
    //=============================================================================
    Scalar ThermoElasticPlaneConstitutive::thick() const
    {
        return _thick;
    }
    //=============================================================================
    ThermoElasticPlaneConstitutive::Type ThermoElasticPlaneConstitutive::type() const
    {
        return _type;
    }
    //=============================================================================
    la::DenseMatrix ThermoElasticPlaneConstitutive::stress_strain_matrix() const
    {
        la::DenseMatrix D(3, 3);
        if (_type == Type::plane_stress)
        {
            Scalar c = _thick * _prop.E / (1 - _prop.nu * _prop.nu);
            D.insert(0, 0, c * 1.0);
            D.insert(0, 1, c * _prop.nu);
            D.insert(1, 0, c * _prop.nu);
            D.insert(1, 1, c * 1.0);
            D.insert(2, 2, c * (1 - _prop.nu) * 0.5);
        }
        else
        {
            Scalar c = _thick * _prop.E / (1 - 2 * _prop.nu) / (1 + _prop.nu);
            D.insert(0, 0, c * (1.0 - _prop.nu));
            D.insert(0, 1, c * _prop.nu);
            D.insert(1, 0, c * _prop.nu);
            D.insert(1, 1, c * (1.0 - _prop.nu));
            D.insert(2, 2, c * (1 - 2 * _prop.nu) * 0.5);
        }
        return D;
    }
    //=============================================================================
    la::DenseMatrix ThermoElasticPlaneConstitutive::strain_displacement_matrix(const std::vector<Scalar> &dNdX) const
    {
        int n_nodes = dNdX.size() / 3;
        int n_cols = n_nodes * 2;
        int n_rows = 3;
        la::DenseMatrix B(n_rows, n_cols);
        for (int i = 0; i < n_nodes; i++)
        {
            // exx
            B.insert(0, i * 2 + 0, dNdX[i * 3 + 0]);

            // eyy
            B.insert(1, i * 2 + 1, dNdX[i * 3 + 1]);

            // exy
            B.insert(2, i * 2 + 0, dNdX[i * 3 + 1]);
            B.insert(2, i * 2 + 1, dNdX[i * 3 + 0]);
        }
        return B;
    }
    //=============================================================================
    Scalar ThermoElasticPlaneConstitutive::eval_vm_stress(const std::vector<Scalar> &dNdx, const std::vector<Scalar> &u) const
    {
        return 0;
    }
    //=============================================================================
    la::DenseMatrix ThermoElasticPlaneConstitutive::eval_thermal_strain(const std::vector<Scalar> &N, const std::vector<Scalar> &dT) const
    {
        Scalar alpha = _prop.alpha;
        Scalar coeff = 1.0;
        if (_type == Type::plane_strain)
        {
            coeff = 1 + _prop.nu;
        }
        la::DenseMatrix eps(3, 1);
        for (std::size_t i = 0; i < N.size(); i++)
        {
            eps.add(0, 0, coeff * alpha * N[i] * dT[i]);
            eps.add(1, 0, coeff * alpha * N[i] * dT[i]);
        }
        return eps;
    }
}