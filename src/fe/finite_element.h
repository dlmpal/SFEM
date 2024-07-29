#pragma once

#include "basis/basis.h"
#include "../la/dense_matrix.h"
#include <array>
#include <memory>

// Forward declaration
namespace sfem::fe
{
    class Function;
}

namespace sfem::fe
{
    /// @brief Element coordinate transform  data
    struct FEData
    {
        /// @brief Quadrature weight
        Scalar qwt = 0;

        /// @brief Quadrature point
        std::array<Scalar, 3> qpt = {0};

        /// @brief Natural-to-physical Jacobian determinant
        Scalar detJ = 0;

        /// @brief Natural-to-physical Jacobian (direct)
        std::array<Scalar, 3 * 3> dXdxi = {0};

        /// @brief Physical-to-natural jacobian (inverse)
        std::array<Scalar, 3 * 3> dxidX = {0};

        /// @brief Shape function
        std::vector<Scalar> N;

        /// @brief Shape function gradient (natural)
        std::vector<Scalar> dNdxi;

        /// @brief Shape function gradient (physical)
        std::vector<Scalar> dNdX;
    };

    /// @brief Matrix types that can be
    /// evaluated by a FiniteElement
    /// @note Not all FEs support all types
    enum class FEMatrixType
    {
        mass = 0,
        damping = 1,
        stiffness = 2,
        jacobian = 3,
    };

    /// @brief Vector types that can be
    /// evaluated by a FiniteElement
    /// @note Not all FEs support all types
    enum class FEVectorType
    {
        load = 0,
        residual = 1
    };

    /// @brief Base FiniteElement class
    class FiniteElement
    {
    public:
        /// @brief Create a FiniteElement
        /// @param name Element name
        /// @param n_vars Number of variables per node
        /// @param physical_dim Element physical dimension
        /// @param cell Corresponding mesh Cell
        FiniteElement(const std::string &name,
                      int n_vars,
                      int physical_dim,
                      mesh::Cell cell);

        /// @brief Virtual destructor
        virtual ~FiniteElement() = default;

        /// @brief Get the element's name
        std::string name() const;

        /// @brief Get the number of variables per node for the element
        int n_vars() const;

        /// @brief Get the element's physical dimension
        int physical_dim() const;

        /// @brief Get the number of DoF for the element
        int n_dof() const;

        /// @brief Get a reference to the element's underlying Cell
        const mesh::Cell &cell() const;

        /// @brief Get a pointer to the element's Basis
        basis::Basis *basis() const;

        /// @brief Evaluate the element's basis transformation at the n-th quadrature point
        virtual FEData transform_basis(int npt, const std::vector<Scalar> &xpts) const;

        /// @brief Evaluate the element mass matrix
        /// @note  Returns a zero matrix if not overwritten
        /// @param data Basis transformation data
        /// @param xpts Element nodal positions
        /// @param u Field values corresponding to the element nodes
        /// @return Element mass matrix
        virtual la::DenseMatrix evaluate_mass_matrix(const FEData &data,
                                                     const std::vector<Scalar> &xpts,
                                                     const std::vector<Scalar> &u) const
        {
            return la::DenseMatrix(n_dof(), n_dof());
        }

        /// @brief Evaluate the element damping matrix
        /// @note  Returns a zero matrix if not overwritten
        /// @param data Basis transformation data
        /// @param xpts Element nodal positions
        /// @param u Field values corresponding to the element nodes
        /// @return Element damping matrix
        virtual la::DenseMatrix evaluate_damping_matrix(const FEData &data,
                                                        const std::vector<Scalar> &xpts,
                                                        const std::vector<Scalar> &u) const
        {
            return la::DenseMatrix(n_dof(), n_dof());
        }

        /// @brief Evaluate the element stiffness matrix
        /// @note  Returns a zero matrix if not overwritten
        /// @param data Basis transformation data
        /// @param xpts Element nodal positions
        /// @param u Field values corresponding to the element nodes
        /// @return Element stiffness matrix
        virtual la::DenseMatrix evaluate_stiff_matrix(const FEData &data,
                                                      const std::vector<Scalar> &xpts,
                                                      const std::vector<Scalar> &u) const
        {
            return la::DenseMatrix(n_dof(), n_dof());
        }

        /// @brief Evaluate the element load vector
        /// @note  Returns a zero vector if not overwritten
        /// @param data Basis transformation data
        /// @param xpts Element nodal positions
        /// @param u Field values corresponding to the element nodes
        /// @return Element load vector
        /// @todo Change return value to std::vector
        virtual la::DenseMatrix evaluate_load_vector(const FEData &data,
                                                     const std::vector<Scalar> &xpts,
                                                     const std::vector<Scalar> &u) const
        {
            return la::DenseMatrix(n_dof(), 1);
        }

        la::DenseMatrix integrate_fe_matrix(const std::vector<Scalar> &xpts,
                                            const std::vector<Scalar> &u,
                                            FEMatrixType type) const;

        la::DenseMatrix integrate_fe_vector(const std::vector<Scalar> &xpts,
                                            const std::vector<Scalar> &u,
                                            FEVectorType type) const;

        la::DenseMatrix integrate_function(const std::vector<Scalar> &xpts,
                                           const std::vector<Scalar> &u,
                                           const Function &func) const;

        std::tuple<la::DenseMatrix, la::DenseMatrix> project_function(const std::vector<Scalar> &xpts,
                                                                      const std::vector<Scalar> &u,
                                                                      const Function &func) const;

    protected:
        /// @brief Element name
        /// @note The name is used to identify the element,
        /// for example when downcasting
        std::string name_;

        /// @brief Number of variables per node
        int n_vars_;

        /// @brief Element physical dimension
        /// @note Not neccessarily equal to the cell's dimension
        int physical_dim_;

        /// @brief Underlying mesh cell
        mesh::Cell cell_;

        /// @brief Element basis
        std::unique_ptr<basis::Basis> basis_;
    };
}