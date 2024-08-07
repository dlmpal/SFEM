#pragma once

#include "../../finite_element.h"

namespace sfem::fe::thermal
{
    class HeatFlux2D : public FiniteElement
    {
    public:
        HeatFlux2D(mesh::Cell cell, Scalar flux, Scalar thick);

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u,
                                             Scalar time = 0) const override;

    private:
        /// @brief Heat flux value
        Scalar flux_;

        /// @brief Thickness
        Scalar thick_;
    };

    class HeatFlux3D : public FiniteElement
    {
    public:
        HeatFlux3D(mesh::Cell cell, Scalar flux);

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u,
                                             Scalar time = 0) const override;

    private:
        /// @brief Heat flux value
        Scalar flux_;
    };
}