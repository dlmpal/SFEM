#include "thermoelastic.h"

namespace sfem::fe::constitutive
{
    //=============================================================================
    ThermoElasticConstitutive::ThermoElasticConstitutive(ThermoMechanicalProperties &prop)
        : _prop(prop)
    {
    }
    //=============================================================================
    ThermoMechanicalProperties &ThermoElasticConstitutive::prop() const
    {
        return _prop;
    }
    //=============================================================================
    la::DenseMatrix ThermoElasticConstitutive::eval_strain(const std::vector<Scalar> &dNdX, const std::vector<Scalar> &u) const
    {
        auto B = strain_displacement_matrix(dNdX);
        auto U = la::DenseMatrix(u.size(), 1, u);
        return B * U;
    }
    //=============================================================================
    la::DenseMatrix ThermoElasticConstitutive::eval_stress(const std::vector<Scalar> &dNdX, const std::vector<Scalar> &u) const
    {
        auto D = stress_strain_matrix();
        auto B = strain_displacement_matrix(dNdX);
        auto U = la::DenseMatrix(u.size(), 1, u);
        return D * B * U;
    }
    //=============================================================================
    la::DenseMatrix ThermoElasticConstitutive::eval_thermal_stress(const std::vector<Scalar> &N, const std::vector<Scalar> &dNdX, const std::vector<Scalar> &dT) const
    {
        auto D = stress_strain_matrix();
        auto B = strain_displacement_matrix(dNdX);
        auto eps = eval_thermal_strain(N, dT);
        return B.T() * D * eps;
    }
}