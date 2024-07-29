#pragma once

#include "function.h"

namespace sfem::fe::function
{
    class Stress2D : public Function
    {
    public:
        Stress2D();

        la::DenseMatrix operator()(const FiniteElement &ele,
                                   const FEData &data,
                                   const std::vector<Scalar> &xpts,
                                   const std::vector<Scalar> &u) const override;
    };

    class VonMisesStress2D : public Function
    {
    public:
        VonMisesStress2D();

        la::DenseMatrix operator()(const FiniteElement &ele,
                                   const FEData &data,
                                   const std::vector<Scalar> &xpts,
                                   const std::vector<Scalar> &u) const override;
    };

    class Stress3D : public Function
    {
    public:
        Stress3D();

        la::DenseMatrix operator()(const FiniteElement &ele,
                                   const FEData &data,
                                   const std::vector<Scalar> &xpts,
                                   const std::vector<Scalar> &u) const override;
    };

    class VonMisesStress3D : public Function
    {
    public:
        VonMisesStress3D();

        la::DenseMatrix operator()(const FiniteElement &ele,
                                   const FEData &data,
                                   const std::vector<Scalar> &xpts,
                                   const std::vector<Scalar> &u) const override;
    };
}