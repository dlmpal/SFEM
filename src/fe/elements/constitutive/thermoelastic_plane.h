#pragma once

#include "thermoelastic.h"

namespace sfem::fe::constitutive
{
    class ThermoElasticPlaneConstitutive : public ThermoElasticConstitutive
    {
    public:
        enum class Type
        {
            plane_stress = 0,
            plane_strain = 1
        };

        /// @brief Create a ThermoElasticPlaneConstitutive
        /// @param prop The corresponding ThermoMechanicalProperties
        /// @param thick Thickness
        /// @param type Plane constitutive type. Either plane stress or plane strain
        ThermoElasticPlaneConstitutive(ThermoMechanicalProperties &prop, Scalar thick, Type = Type::plane_stress);

        /// @brief Get the thickness
        Scalar thick() const;

        /// @brief Get the plane constitutive type
        Type type() const;

        /// @brief Evaluate the stress-strain matrix
        la::DenseMatrix stress_strain_matrix() const override;

        /// @brief Evaluate the strain-displacement matrix
        la::DenseMatrix strain_displacement_matrix(const std::vector<Scalar> &dNdX) const override;

        /// @brief Evaluate the von Mises stress, given the displacement
        Scalar eval_vm_stress(const std::vector<Scalar> &dNdx, const std::vector<Scalar> &u) const;

        /// @brief Evaluate the thermal strain, given the temperature difference
        la::DenseMatrix eval_thermal_strain(const std::vector<Scalar> &N, const std::vector<Scalar> &dT) const override;

    private:
        /// @brief Thickness
        Scalar _thick;

        /// @brief Type
        Type _type;
    };
}