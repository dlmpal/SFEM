#pragma once

#include "../../finite_element.h"
#include "../../../mesh/field.h"

namespace sfem::fe::solid
{
    class PressureLoad2D : public FiniteElement
    {
    public:
        PressureLoad2D(mesh::Cell cell, mesh::Field &P, Scalar thick);

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u) const override;

    private:
        /// @brief Pressure field
        mesh::Field &P_;

        /// @brief Thickness
        Scalar thick_;
    };

    class PressureLoad3D : public FiniteElement
    {
    public:
        PressureLoad3D(mesh::Cell cell, mesh::Field &P);

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u) const override;

    private:
        /// @brief Pressure field
        mesh::Field &P_;
    };
}