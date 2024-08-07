#include "heat_flux.h"

namespace sfem::fe::thermal
{
    //=============================================================================
    HeatFlux2D::HeatFlux2D(mesh::Cell cell, Scalar flux, Scalar thick)
        : FiniteElement("HeatFlux2D", 1, 2, cell), flux_(flux), thick_(thick)
    {
    }
    //=============================================================================
    la::DenseMatrix HeatFlux2D::evaluate_load_vector(const FEData &data,
                                                     const std::vector<Scalar> &xpts,
                                                     const std::vector<Scalar> &u,
                                                     Scalar time) const
    {
        int n_nodes = basis_->n_nodes();
        la::DenseMatrix Fe(n_nodes, 1);
        for (int i = 0; i < n_nodes; i++)
        {
            Fe.add(i, 0, thick_ * flux_ * data.N[i]);
        }
        return Fe;
    }
    //=============================================================================
    HeatFlux3D::HeatFlux3D(mesh::Cell cell, Scalar flux)
        : FiniteElement("HeatFlux3D", 1, 3, cell), flux_(flux)
    {
    }
    //=============================================================================
    la::DenseMatrix HeatFlux3D::evaluate_load_vector(const FEData &data,
                                                     const std::vector<Scalar> &xpts,
                                                     const std::vector<Scalar> &u,
                                                     Scalar time) const
    {
        int n_nodes = basis_->n_nodes();
        la::DenseMatrix Fe(n_nodes, 1);
        for (int i = 0; i < n_nodes; i++)
        {
            Fe.add(i, 0, flux_ * data.N[i]);
        }
        return Fe;
    }
}