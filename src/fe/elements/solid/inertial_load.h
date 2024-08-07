#pragma once

#include "../../finite_element.h"
#include "../constitutive/thermoelastic_plane.h"
#include "../constitutive/thermoelastic_solid.h"

namespace sfem::fe::solid
{
    class InertialLoad2D : public FiniteElement
    {
    public:
        InertialLoad2D(const mesh::Cell cell,
                       constitutive::ThermoElasticPlaneConstitutive &constitutive,
                       const std::array<Scalar, 2> &g);

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u,
                                             Scalar time = 0) const override;

    private:
        constitutive::ThermoElasticPlaneConstitutive &constitutive_;

        /// @brief The value of the inertial acceleration
        std::array<Scalar, 2> g_;
    };

    class InertialLoad3D : public FiniteElement
    {
    public:
        InertialLoad3D(const mesh::Cell cell,
                       constitutive::ThermoElasticSolidConstitutive &constitutive,
                       const std::array<Scalar, 3> &g);

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u,
                                             Scalar time = 0) const override;

    private:
        constitutive::ThermoElasticSolidConstitutive &constitutive_;

        /// @brief The value of the inertial acceleration
        std::array<Scalar, 3> g_;
    };
}