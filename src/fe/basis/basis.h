#pragma once

#include "../../mesh/cell.h"

namespace sfem::fe::basis
{
    /// @brief Base Basis class
    class Basis
    {
    public:
        /// @brief Destructor
        virtual ~Basis() = default;

        /// @brief Get the basis dimension
        virtual int dim() const = 0;

        /// @brief Get the number of nodes
        virtual int n_nodes() const = 0;

        /// @brief Get the number of quadrature points
        virtual int n_qpts() const = 0;

        /// @brief Get the quadrature weight
        virtual Scalar qwt(int npt) const = 0;

        /// @brief Get the quadrature point coordiniates
        virtual void qpt(int npt, Scalar pt[]) const = 0;

        /// @brief Evaluate shape function @ pt
        virtual void eval_shape(const Scalar pt[], Scalar N[]) const = 0;

        /// @brief Evaluate shape function gradient @ pt
        virtual void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const = 0;
    };

    /// @brief Point
    class PointBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief Linear 2-node line segment
    class LinearLineBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief  Quadratic 3-node line segment
    class QuadraticLineBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief  Cubic 4-node line segment
    class CubicLineBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief Linear 3-node triangle
    class LinearTriangleBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief  Quadratic 6-node triangle
    class QuadraticTriangleBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief  Cubic 10-node triangle
    /// @todo fix
    class CubicTriangleBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief Linear 4-node quadrilateral
    class LinearQuadBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief Quadratic 9-node quadrilateral
    class QuadraticQuadBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief Cubic 16-node quadrilateral
    class CubicQuadBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief Linear 4-node tetrahedron
    class LinearTetrahedralBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief  Quadratic 10-node tetrahedron
    class QuadraticTehtrahedralBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief  Cubic 20-node tetrahedron
    /// @todo fix
    class CubicTetrahedralBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief  Linear 8-node hexahedron
    class LinearHexahedralBasis : public Basis
    {
    public:
        int dim() const override;
        int n_nodes() const override;
        int n_qpts() const override;
        Scalar qwt(int npt) const override;
        void qpt(int npt, Scalar pt[]) const override;
        void eval_shape(const Scalar pt[], Scalar N[]) const override;
        void eval_shape_grad(const Scalar pt[], Scalar dNdxi[]) const override;
    };

    /// @brief
    Basis *CreateBasis(const mesh::Cell &cell);
}