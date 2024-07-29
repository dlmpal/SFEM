#pragma once

#include "thermoelastic.h"

namespace sfem::fe::constitutive
{
    class ThermoElasticSolidConstitutive : public ThermoElasticConstitutive
    {
    public:
        /// @brief Create a ThermoElasticSolidConstitutive
        /// @param prop The corresponding ThermoMechanicalProperties
        ThermoElasticSolidConstitutive(ThermoMechanicalProperties &prop);

        /// @brief Evaluate the stress-strain matrix
        la::DenseMatrix stress_strain_matrix() const override;

        /// @brief Evaluate the strain-displacement matrix
        la::DenseMatrix strain_displacement_matrix(const std::vector<Scalar> &dNdX) const override;

        /// @brief Evaluate the von Mises stress, given the displacement
        Scalar eval_vm_stress(const std::vector<Scalar> &dNdx, const std::vector<Scalar> &u) const;

        /// @brief Evaluate the thermal strain, given the temperature difference
        la::DenseMatrix eval_thermal_strain(const std::vector<Scalar> &N, const std::vector<Scalar> &dT) const override;
    };
}