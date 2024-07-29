#pragma once

#include "../../finite_element.h"
#include "../../../mesh/field.h"

namespace sfem::fe::thermal
{
    class HeatConvection2D : public FiniteElement
    {
    public:
        HeatConvection2D(mesh::Cell cell, Scalar htc, Scalar T_bulk, Scalar thick);

        la::DenseMatrix evaluate_stiff_matrix(const FEData &data,
                                              const std::vector<Scalar> &xpts,
                                              const std::vector<Scalar> &u) const override;

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u) const override;

    private:
        /// @brief Heat transfer coefficient
        Scalar htc_;

        /// @brief Bulk temperature
        Scalar T_bulk_;

        /// @brief Thickness
        Scalar thick_;
    };

    class HeatConvection3D : public FiniteElement
    {
    public:
        HeatConvection3D(mesh::Cell cell, Scalar htc, Scalar T_bulk);

        la::DenseMatrix evaluate_stiff_matrix(const FEData &data,
                                              const std::vector<Scalar> &xpts,
                                              const std::vector<Scalar> &u) const override;

        la::DenseMatrix evaluate_load_vector(const FEData &data,
                                             const std::vector<Scalar> &xpts,
                                             const std::vector<Scalar> &u) const override;

    private:
        /// @brief Heat transfer coefficient
        Scalar htc_;

        /// @brief Bulk temperature
        Scalar T_bulk_;
    };
}