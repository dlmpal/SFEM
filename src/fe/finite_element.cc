#include "finite_element.h"
#include "functions/function.h"
#include "../common/error.h"
#include "../common/math.h"

namespace sfem::fe
{
    //=============================================================================
    FiniteElement::FiniteElement(const std::string &name,
                                 int n_vars,
                                 int physical_dim,
                                 mesh::Cell cell)
        : name_(name),
          n_vars_(n_vars),
          physical_dim_(physical_dim),
          cell_(cell)
    {
        basis_ = std::unique_ptr<basis::Basis>(basis::CreateBasis(cell));
    }
    //=============================================================================
    std::string FiniteElement::name() const
    {
        return name_;
    }
    //=============================================================================
    int FiniteElement::n_vars() const
    {
        return n_vars_;
    }
    //=============================================================================
    int FiniteElement::physical_dim() const
    {
        return physical_dim_;
    }
    //=============================================================================
    int FiniteElement::n_dof() const
    {
        return basis_->n_nodes() * n_vars_;
    }
    //=============================================================================
    const mesh::Cell &FiniteElement::cell() const
    {
        return cell_;
    }
    //=============================================================================
    basis::Basis *FiniteElement::basis() const
    {
        return basis_.get();
    }
    //=============================================================================
    FEData FiniteElement::transform_basis(int npt, const std::vector<Scalar> &xpts) const
    {
        // Transform data
        FEData data;

        // Allocate memory for the shape function and shape function gradient vectors
        data.N.resize(basis_->n_nodes(), 0);
        data.dNdxi.resize(basis_->n_nodes() * 3, 0);
        data.dNdX.resize(basis_->n_nodes() * 3, 0);

        // Return early for point elements
        if (physical_dim() == 0)
        {
            data.qwt = 1.0;
            data.detJ = 1.0;
            return data;
        }

        // Get the quadrature weight and point
        data.qwt = basis_->qwt(npt);
        basis_->qpt(npt, data.qpt.data());

        // Evaluate the shape function and its gradient w.r.t natural coordinates
        basis_->eval_shape(data.qpt.data(), data.N.data());
        basis_->eval_shape_grad(data.qpt.data(), data.dNdxi.data());

        // Evaluate the natural to physical Jacobian
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < basis_->n_nodes(); k++)
                {
                    data.dXdxi[i * 3 + j] += data.dNdxi[k * 3 + j] * xpts[k * 3 + i];
                }
            }
        }

        // Evaluate the physical to natural, or inverse, Jacobian
        // and the determinant
        if (physical_dim() == basis_->dim())
        {
            data.detJ = math::inv(physical_dim(), data.dXdxi.data(), data.dxidX.data());
        }
        else
        {
            data.detJ = math::pinv(physical_dim(), basis_->dim(), data.dXdxi.data(), data.dxidX.data());
        }

        // Check for non-positive jacobian
        if (data.detJ <= 0)
        {
            error::negative_jacobian_error(cell_.idx(), __FILE__, __LINE__);
        }

        // Evaluate the shape function gradient w.r.t physical coordinates
        math::matmult(basis_->n_nodes(), 3, 3, data.dNdxi.data(), data.dxidX.data(), data.dNdX.data());

        return data;
    }
    //=============================================================================
    la::DenseMatrix FiniteElement::integrate_fe_matrix(const std::vector<Scalar> &xpts,
                                                       const std::vector<Scalar> &u,
                                                       FEMatrixType type) const
    {
        la::DenseMatrix M(n_dof(), n_dof());

        for (int npt = 0; npt < basis_->n_qpts(); npt++)
        {
            auto data = transform_basis(npt, xpts);
            Scalar qwt_detJ = data.qwt * data.detJ;
            switch (type)
            {
            case FEMatrixType::mass:
                M += evaluate_mass_matrix(data, xpts, u) * qwt_detJ;
                break;
            case FEMatrixType::damping:
                M += evaluate_damping_matrix(data, xpts, u) * qwt_detJ;
                break;
            case FEMatrixType::stiffness:
                M += evaluate_stiff_matrix(data, xpts, u) * qwt_detJ;
                break;
            default:
                break;
            }
        }

        return M;
    }
    //=============================================================================
    la::DenseMatrix FiniteElement::integrate_fe_vector(const std::vector<Scalar> &xpts,
                                                       const std::vector<Scalar> &u,
                                                       FEVectorType type) const
    {
        la::DenseMatrix F(n_dof(), 1);

        for (int npt = 0; npt < basis_->n_qpts(); npt++)
        {
            auto data = transform_basis(npt, xpts);
            Scalar qwt_detJ = data.qwt * data.detJ;
            switch (type)
            {
            case FEVectorType::load:
                F += evaluate_load_vector(data, xpts, u) * qwt_detJ;
                break;
            default:
                break;
            }
        }

        return F;
    }
    //=============================================================================
    la::DenseMatrix FiniteElement::integrate_function(const std::vector<Scalar> &xpts,
                                                      const std::vector<Scalar> &u,
                                                      const Function &func) const
    {
        la::DenseMatrix F(func.size(), 1);
        for (int npt = 0; npt < basis_->n_qpts(); npt++)
        {
            auto data = transform_basis(npt, xpts);
            F += func(*this, data, xpts, u) * data.detJ * data.qwt;
        }
        return F;
    }
    //=============================================================================
    std::tuple<la::DenseMatrix, la::DenseMatrix>
    FiniteElement::project_function(const std::vector<Scalar> &xpts,
                                    const std::vector<Scalar> &u,
                                    const Function &func) const
    {
        la::DenseMatrix M(cell_.n_nodes(), cell_.n_nodes());
        la::DenseMatrix F(cell_.n_nodes(), func.size());
        for (int npt = 0; npt < basis_->n_qpts(); npt++)
        {
            auto data = transform_basis(npt, xpts);
            auto values = func(*this, data, xpts, u);

            for (int i = 0; i < cell_.n_nodes(); i++)
            {
                for (int j = 0; j < cell_.n_nodes(); j++)
                {
                    M.add(i, j, data.N[i] * data.N[j] * data.detJ * data.qwt);
                }
                for (int j = 0; j < func.size(); j++)
                {
                    F.add(i, j, data.N[i] * values.at(j, 0) * data.detJ * data.qwt);
                }
            }
        }

        return std::make_tuple(M, F);
    }
}