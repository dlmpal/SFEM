#include "heat_conduction.h"

namespace sfem::fe::thermal
{
    //=============================================================================
    HeatConduction2D::HeatConduction2D(const mesh::Cell cell,
                                       constitutive::ThermoElasticPlaneConstitutive &constitutive)
        : FiniteElement("HeatConduction2D", 1, 2, cell),
          constitutive_(constitutive)
    {
    }
    //=============================================================================
    constitutive::ThermoElasticPlaneConstitutive &HeatConduction2D::constitutive() const
    {
        return constitutive_;
    }
    //=============================================================================
    void HeatConduction2D::add_heat_load(const std::shared_ptr<FiniteElement> &load)
    {
        loads_.push_back(load);
    }
    //=============================================================================
    la::DenseMatrix HeatConduction2D::evaluate_mass_matrix(const FEData &data,
                                                           const std::vector<Scalar> &xpt,
                                                           const std::vector<Scalar> &u,
                                                           Scalar time) const
    {
        int n_nodes = basis_->n_nodes();
        Scalar rho = constitutive_.prop().rho;
        Scalar cp = constitutive_.prop().cp;
        Scalar thick = constitutive_.thick();
        la::DenseMatrix Me(n_nodes, n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            for (int j = 0; j < n_nodes; j++)
            {
                Me.add(i, j, rho * cp * thick * data.N[i] * data.N[j]);
            }
        }
        return Me;
    }
    //=============================================================================
    la::DenseMatrix HeatConduction2D::evaluate_stiff_matrix(const FEData &data,
                                                            const std::vector<Scalar> &xpts,
                                                            const std::vector<Scalar> &u,
                                                            Scalar time) const
    {
        int n_nodes = basis_->n_nodes();
        Scalar kappa = constitutive_.prop().kappa;
        Scalar thick = constitutive_.thick();
        la::DenseMatrix Ke(n_nodes, n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            for (int j = 0; j < n_nodes; j++)
            {
                Ke.add(i, j, thick * kappa * data.dNdX[i * 3 + 0] * data.dNdX[j * 3 + 0]);
                Ke.add(i, j, thick * kappa * data.dNdX[i * 3 + 1] * data.dNdX[j * 3 + 1]);
            }
        }
        return Ke;
    }
    //=============================================================================
    la::DenseMatrix HeatConduction2D::evaluate_load_vector(const FEData &data,
                                                           const std::vector<Scalar> &xpts,
                                                           const std::vector<Scalar> &u,
                                                           Scalar time) const
    {
        la::DenseMatrix Fe(basis_->n_nodes(), 1);
        for (std::size_t i = 0; i < loads_.size(); i++)
        {
            Fe += loads_[i]->evaluate_load_vector(data, xpts, u);
        }
        return Fe;
    }
    //=============================================================================
    HeatConduction3D::HeatConduction3D(const mesh::Cell cell,
                                       constitutive::ThermoElasticSolidConstitutive &constitutive)
        : FiniteElement("HeatConduction3D", 1, 3, cell),
          constitutive_(constitutive)
    {
    }
    //=============================================================================
    constitutive::ThermoElasticSolidConstitutive &HeatConduction3D::constitutive() const
    {
        return constitutive_;
    }
    //=============================================================================
    void HeatConduction3D::add_heat_load(const std::shared_ptr<FiniteElement> &load)
    {
        loads_.push_back(std::shared_ptr<FiniteElement>(load));
    }
    //=============================================================================
    la::DenseMatrix HeatConduction3D::evaluate_mass_matrix(const FEData &data,
                                                           const std::vector<Scalar> &xpts,
                                                           const std::vector<Scalar> &u,
                                                           Scalar time) const
    {
        int n_nodes = basis_->n_nodes();
        Scalar rho = constitutive_.prop().rho;
        Scalar cp = constitutive_.prop().cp;
        la::DenseMatrix Me(n_nodes, n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            for (int j = 0; j < n_nodes; j++)
            {
                Me.add(i, j, rho * cp * data.N[i] * data.N[j]);
            }
        }
        return Me;
    }
    //=============================================================================
    la::DenseMatrix HeatConduction3D::evaluate_stiff_matrix(const FEData &data,
                                                            const std::vector<Scalar> &xpts,
                                                            const std::vector<Scalar> &u,
                                                            Scalar time) const
    {
        int n_nodes = basis_->n_nodes();
        Scalar kappa = constitutive_.prop().kappa;
        la::DenseMatrix Ke(n_nodes, n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            for (int j = 0; j < n_nodes; j++)
            {
                Ke.add(i, j, kappa * data.dNdX[i * 3 + 0] * data.dNdX[j * 3 + 0]);
                Ke.add(i, j, kappa * data.dNdX[i * 3 + 1] * data.dNdX[j * 3 + 1]);
                Ke.add(i, j, kappa * data.dNdX[i * 3 + 2] * data.dNdX[j * 3 + 2]);
            }
        }
        return Ke;
    }
    //=============================================================================
    la::DenseMatrix HeatConduction3D::evaluate_load_vector(const FEData &data,
                                                           const std::vector<Scalar> &xpts,
                                                           const std::vector<Scalar> &u,
                                                           Scalar time) const
    {
        la::DenseMatrix Fe(basis_->n_nodes(), 1);
        for (std::size_t i = 0; i < loads_.size(); i++)
        {
            Fe += loads_[i]->evaluate_load_vector(data, xpts, u);
        }
        return Fe;
    }
}