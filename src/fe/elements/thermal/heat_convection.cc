#include "heat_convection.h"

namespace sfem::fe::thermal
{
    //=============================================================================
    HeatConvection2D::HeatConvection2D(mesh::Cell cell, Scalar htc, Scalar T_bulk, Scalar thick)
        : FiniteElement("HeatFlux2D", 1, 2, cell), htc_(htc), T_bulk_(T_bulk), thick_(thick)
    {
    }
    //=============================================================================
    la::DenseMatrix HeatConvection2D::evaluate_stiff_matrix(const FEData &data,
                                                            const std::vector<Scalar> &xpts,
                                                            const std::vector<Scalar> &u,
                                                            Scalar time = 0) const
    {
        int n_nodes = basis_->n_nodes();
        la::DenseMatrix Ke(n_nodes, n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            for (int j = 0; j < n_nodes; j++)
            {
                Ke.add(i, j, thick_ * htc_ * data.N[i] * data.N[j]);
            }
        }
        return Ke;
    }
    //=============================================================================
    la::DenseMatrix HeatConvection2D::evaluate_load_vector(const FEData &data,
                                                           const std::vector<Scalar> &xpts,
                                                           const std::vector<Scalar> &u,
                                                           Scalar time = 0) const
    {
        int n_nodes = basis_->n_nodes();
        la::DenseMatrix Fe(n_nodes, 1);
        for (int i = 0; i < n_nodes; i++)
        {
            Fe.add(i, 0, thick_ * htc_ * T_bulk_ * data.N[i]);
        }
        return Fe;
    }
    //=============================================================================
    HeatConvection3D::HeatConvection3D(mesh::Cell cell, Scalar htc, Scalar T_bulk)
        : FiniteElement("HeatFlux2D", 1, 3, cell), htc_(htc), T_bulk_(T_bulk)
    {
    }
    //=============================================================================
    la::DenseMatrix HeatConvection3D::evaluate_stiff_matrix(const FEData &data,
                                                            const std::vector<Scalar> &xpts,
                                                            const std::vector<Scalar> &u,
                                                            Scalar time = 0) const
    {
        int n_nodes = basis_->n_nodes();
        la::DenseMatrix Ke(n_nodes, n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            for (int j = 0; j < n_nodes; j++)
            {
                Ke.add(i, j, htc_ * data.N[i] * data.N[j]);
            }
        }
        return Ke;
    }
    //=============================================================================
    la::DenseMatrix HeatConvection3D::evaluate_load_vector(const FEData &data,
                                                           const std::vector<Scalar> &xpts,
                                                           const std::vector<Scalar> &u,
                                                           Scalar time = 0) const
    {
        int n_nodes = basis_->n_nodes();
        la::DenseMatrix Fe(n_nodes, 1);
        for (int i = 0; i < n_nodes; i++)
        {
            Fe.add(i, 0, htc_ * T_bulk_ * data.N[i]);
        }
        return Fe;
    }
}