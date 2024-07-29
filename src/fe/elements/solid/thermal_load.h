#pragma once

#include "../../finite_element.h"
#include "../constitutive/thermoelastic_plane.h"
#include "../constitutive/thermoelastic_solid.h"
#include "../../../mesh/field.h"

namespace sfem::fe::solid
{
    class ThermalLoad2D : fe::FiniteElement
    {
    public:
        ThermalLoad2D(mesh::Cell cell,
                      constitutive::ThermoElasticPlaneConstitutive &constitutive,
                      mesh::Field &T,
                      Scalar T0);

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u) const override;

    private:
        /// @brief Constitutive
        constitutive::ThermoElasticPlaneConstitutive &constitutive_;

        /// @brief Temperature field
        const mesh::Field &T_;

        /// @brief Initial temperature
        Scalar T0_;
    };
}