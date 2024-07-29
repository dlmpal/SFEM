#pragma once

#include "function.h"
#include "../../mesh/field.h"

namespace sfem::fe::function
{
    class FieldGradient : public Function
    {
    public:
        FieldGradient(const mesh::Field &field);

        la::DenseMatrix operator()(const FiniteElement &elem,
                                   const FEData &data,
                                   const std::vector<Scalar> &xpts,
                                   const std::vector<Scalar> &u) const override;
    };
}