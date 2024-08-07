#pragma once

#include "function.h"

namespace sfem::fe::function
{
    class StructuralMass2D : public Function
    {
    public:
        StructuralMass2D();

        la::DenseMatrix operator()(const FiniteElement &ele,
                                   const FEData &data,
                                   const std::vector<Scalar> &xpts,
                                   const std::vector<Scalar> &u,
                                   Scalar time = 0) const override;
    };

    class StructuralMass3D : public Function
    {
    public:
        StructuralMass3D();

        la::DenseMatrix operator()(const FiniteElement &ele,
                                   const FEData &data,
                                   const std::vector<Scalar> &xpts,
                                   const std::vector<Scalar> &u,
                                   Scalar time = 0) const override;
    };
}