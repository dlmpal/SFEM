#pragma once

#include "../../finite_element.h"
#include "../constitutive/thermoelastic_plane.h"
#include "../constitutive/thermoelastic_solid.h"
#include <vector>

namespace sfem::fe::thermal
{
    class HeatConduction2D : public FiniteElement
    {
    public:
        HeatConduction2D(const mesh::Cell cell, constitutive::ThermoElasticPlaneConstitutive &constitutive);

        constitutive::ThermoElasticPlaneConstitutive &constitutive() const;

        void add_heat_load(const std::shared_ptr<FiniteElement> &load);

        la::DenseMatrix evaluate_mass_matrix(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u) const override;

        la::DenseMatrix evaluate_stiff_matrix(const FEData &data,
                                              const std::vector<Scalar> &xpts,
                                              const std::vector<Scalar> &u) const override;

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u) const override;

    private:
        /// @brief Constitutive
        constitutive::ThermoElasticPlaneConstitutive &constitutive_;

        /// @brief List of body loads
        std::vector<std::shared_ptr<FiniteElement>> loads_;
    };

    class HeatConduction3D : public FiniteElement
    {
    public:
        HeatConduction3D(const mesh::Cell cell, constitutive::ThermoElasticSolidConstitutive &constitutive);

        constitutive::ThermoElasticSolidConstitutive &constitutive() const;

        void add_heat_load(const std::shared_ptr<FiniteElement> &load);

        la::DenseMatrix evaluate_mass_matrix(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u) const override;

        la::DenseMatrix evaluate_stiff_matrix(const FEData &data,
                                              const std::vector<Scalar> &xpts,
                                              const std::vector<Scalar> &u) const override;

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u) const override;

    private:
        /// @brief Constitutive
        constitutive::ThermoElasticSolidConstitutive &constitutive_;

        /// @brief List of body loads
        std::vector<std::shared_ptr<FiniteElement>> loads_;
    };
}