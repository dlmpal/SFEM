#pragma once

#include "../../finite_element.h"
#include "../constitutive/thermoelastic_plane.h"
#include "../constitutive/thermoelastic_solid.h"
#include "../../../mesh/field.h"

namespace sfem::fe::solid
{
    class LinearElasticity2D : public FiniteElement
    {
    public:
        LinearElasticity2D(mesh::Cell cell, constitutive::ThermoElasticPlaneConstitutive &constitutive);

        constitutive::ThermoElasticPlaneConstitutive &constitutive() const;

        void add_inertial_load(const std::array<Scalar, 2> &g);

        void add_thermal_load(mesh::Field &T, Scalar T0);

        la::DenseMatrix evaluate_mass_matrix(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u,
                                             Scalar time = 0) const override;

        la::DenseMatrix evaluate_stiff_matrix(const FEData &data,
                                              const std::vector<Scalar> &xpts,
                                              const std::vector<Scalar> &u,
                                              Scalar time = 0) const override;

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u,
                                             Scalar time = 0) const override;

    private:
        /// @brief Constitutive
        constitutive::ThermoElasticPlaneConstitutive &constitutive_;

        /// @brief List of body loads
        std::vector<std::shared_ptr<FiniteElement>> loads_;
    };

    class LinearElasticity3D : public FiniteElement
    {
    public:
        LinearElasticity3D(mesh::Cell cell, constitutive::ThermoElasticSolidConstitutive &constitutive);

        constitutive::ThermoElasticSolidConstitutive &constitutive() const;

        void add_inertial_load(const std::array<Scalar, 3> &g);

        la::DenseMatrix evaluate_mass_matrix(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u,
                                             Scalar time = 0) const override;

        la::DenseMatrix evaluate_stiff_matrix(const FEData &data,
                                              const std::vector<Scalar> &xpts,
                                              const std::vector<Scalar> &u,
                                              Scalar time = 0) const override;

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u,
                                             Scalar time = 0) const override;

    private:
        /// @brief Constitutive
        constitutive::ThermoElasticSolidConstitutive &constitutive_;

        /// @brief List of body loads
        std::vector<std::shared_ptr<FiniteElement>> loads_;
    };
}