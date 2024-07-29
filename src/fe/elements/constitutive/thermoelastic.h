#pragma once

#include "../../../la/dense_matrix.h"

namespace sfem::fe::constitutive
{
    struct ThermoMechanicalProperties
    {
        /// @brief Young's modulus
        Scalar E;

        /// @brief Poisson's ratio
        Scalar nu;

        /// @brief Thermal conductivity
        Scalar kappa;

        /// @brief Density
        Scalar rho;

        /// @brief Specific heat
        Scalar cp;

        /// @brief Heat expansion coefficient
        Scalar alpha;
    };

    /// @brief Base ThermoElasticConstitutive class
    class ThermoElasticConstitutive
    {
    public:
        /// @brief Create a ThermoElasticConstitutive
        /// @param prop The corresponding ThermoMechanicalProperties
        ThermoElasticConstitutive(ThermoMechanicalProperties &prop);

        /// @brief Virtual destructor
        virtual ~ThermoElasticConstitutive() = default;

        /// @brief Get a reference to the proprties
        ThermoMechanicalProperties &prop() const;

        /// @brief Evaluate the stress-strain matrix
        virtual la::DenseMatrix stress_strain_matrix() const = 0;

        /// @brief Evaluate the strain-displacement matrix
        /// @param dNdX Shape function gradient
        virtual la::DenseMatrix strain_displacement_matrix(const std::vector<Scalar> &dNdX) const = 0;

        /// @brief Evaluate the strain, given the displacement
        /// @param dNdX Shape function gradient
        /// @param u Displacement
        /// @return Strain
        la::DenseMatrix eval_strain(const std::vector<Scalar> &dNdX, const std::vector<Scalar> &u) const;

        /// @brief Evaluate the stress, given the displacement
        /// @param dNdX Shape function gradient
        /// @param u Displacement
        /// @return Stress
        la::DenseMatrix eval_stress(const std::vector<Scalar> &dNdX, const std::vector<Scalar> &u) const;

        /// @brief Evaluate the von Mises stress, given the displacement
        /// @param dNdX Shape function gradient
        /// @param u Displacement
        /// @return von Mises stress
        virtual Scalar eval_vm_stress(const std::vector<Scalar> &dNdx, const std::vector<Scalar> &u) const = 0;

        /// @brief Evaluate the thermal strain, given the temperature difference
        /// @param N Shape function
        /// @param dT Temperature difference
        /// @return Thermal strain
        virtual la::DenseMatrix eval_thermal_strain(const std::vector<Scalar> &N, const std::vector<Scalar> &dT) const = 0;

        /// @brief Evaluate the thermal stress, given the temperature difference
        /// @param N Shape function
        /// @param dNdX Shape function gradient
        /// @param dT Temperature difference
        /// @return Thermal stress
        la::DenseMatrix eval_thermal_stress(const std::vector<Scalar> &N, const std::vector<Scalar> &dNdX, const std::vector<Scalar> &dT) const;

    protected:
        ThermoMechanicalProperties &_prop;
    };
}