#pragma once

#include "../../la/dense_matrix.h"

// Forward declaration
namespace sfem::fe
{
    class FEData;
    class FiniteElement;
}

namespace sfem::fe
{
    /// @brief Base function class
    class Function
    {
    public:
        Function(int size)
            : size_(size)
        {
        }

        /// @brief Get the output size
        int size() const
        {
            return size_;
        }

        virtual la::DenseMatrix operator()(const FiniteElement &elem,
                                           const FEData &data,
                                           const std::vector<Scalar> &xpts,
                                           const std::vector<Scalar> &u,
                                           Scalar time = 0) const = 0;

    protected:
        int size_;
    };
}