#include "linear_elasticity.h"
#include "inertial_load.h"
#include "thermal_load.h"

namespace sfem::fe::solid
{
    //=============================================================================
    LinearElasticity2D::LinearElasticity2D(mesh::Cell cell,
                                           constitutive::ThermoElasticPlaneConstitutive &constitutive)
        : FiniteElement("LinearElasticity2D", 2, 2, cell),
          constitutive_(constitutive)
    {
    }
    //=============================================================================
    constitutive::ThermoElasticPlaneConstitutive &LinearElasticity2D::constitutive() const
    {
        return constitutive_;
    }
    //=============================================================================
    void LinearElasticity2D::add_inertial_load(const std::array<Scalar, 2> &g)
    {
        auto load = std::make_shared<InertialLoad2D>(cell_, constitutive_, g);
        loads_.push_back(load);
    }
    //=============================================================================
    void LinearElasticity2D::add_thermal_load(mesh::Field &T, Scalar T0)
    {
        // auto load = std::make_shared<ThermalLoad2D>(cell_, constitutive_, T, T0);
        // loads_.push_back(load);
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity2D::evaluate_mass_matrix(const FEData &data,
                                                             const std::vector<Scalar> &xpts,
                                                             const std::vector<Scalar> &u,
                                                             Scalar time) const
    {
        la::DenseMatrix Me(n_dof(), n_dof());
        Scalar thick = constitutive_.thick();
        Scalar rho = constitutive_.prop().rho;
        for (int i = 0; i < basis_->n_nodes(); i++)
        {
            for (int j = 0; j < basis_->n_nodes(); j++)
            {
                Me.add(i * 2 + 0, j * 2 + 0, thick * rho * data.N[i] * data.N[j]);
                Me.add(i * 2 + 1, j * 2 + 1, thick * rho * data.N[i] * data.N[j]);
            }
        }
        return Me;
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity2D::evaluate_stiff_matrix(const FEData &data,
                                                              const std::vector<Scalar> &xpts,
                                                              const std::vector<Scalar> &u,
                                                              Scalar time) const
    {
        auto B = constitutive_.strain_displacement_matrix(data.dNdX);
        auto D = constitutive_.stress_strain_matrix();
        return B.T() * D * B;
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity2D::evaluate_load_vector(const FEData &data,
                                                             const std::vector<Scalar> &xpts,
                                                             const std::vector<Scalar> &u,
                                                             Scalar time) const
    {
        la::DenseMatrix Fe(n_dof(), 1);
        for (std::size_t i = 0; i < loads_.size(); i++)
        {
            Fe += loads_[i]->evaluate_load_vector(data, xpts, u);
        }
        return Fe;
    }
    //=============================================================================
    LinearElasticity3D::LinearElasticity3D(mesh::Cell cell,
                                           constitutive::ThermoElasticSolidConstitutive &constitutive)
        : FiniteElement("LinearElasticity3D", 3, 3, cell),
          constitutive_(constitutive)
    {
    }
    //=============================================================================
    constitutive::ThermoElasticSolidConstitutive &LinearElasticity3D::constitutive() const
    {
        return constitutive_;
    }
    //=============================================================================
    void LinearElasticity3D::add_inertial_load(const std::array<Scalar, 3> &g)
    {
        auto load = std::make_shared<InertialLoad3D>(cell_, constitutive_, g);
        loads_.push_back(load);
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity3D::evaluate_stiff_matrix(const FEData &data,
                                                              const std::vector<Scalar> &xpts,
                                                              const std::vector<Scalar> &u,
                                                              Scalar time) const
    {
        auto B = constitutive_.strain_displacement_matrix(data.dNdX);
        auto D = constitutive_.stress_strain_matrix();
        return B.T() * D * B;
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity3D::evaluate_mass_matrix(const FEData &data,
                                                             const std::vector<Scalar> &xpts,
                                                             const std::vector<Scalar> &u,
                                                             Scalar time) const
    {
        la::DenseMatrix Me(n_dof(), n_dof());
        Scalar rho = constitutive_.prop().rho;
        for (int i = 0; i < basis_->n_nodes(); i++)
        {
            for (int j = 0; j < basis_->n_nodes(); j++)
            {
                Me.add(i * 3 + 0, j * 3 + 0, rho * data.N[i] * data.N[j]);
                Me.add(i * 3 + 1, j * 3 + 1, rho * data.N[i] * data.N[j]);
                Me.add(i * 3 + 2, j * 3 + 2, rho * data.N[i] * data.N[j]);
            }
        }
        return Me;
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity3D::evaluate_load_vector(const FEData &data,
                                                             const std::vector<Scalar> &xpts,
                                                             const std::vector<Scalar> &u,
                                                             Scalar time) const
    {
        la::DenseMatrix Fe(n_dof(), 1);
        for (std::size_t i = 0; i < loads_.size(); i++)
        {
            Fe += loads_[i]->evaluate_load_vector(data, xpts, u);
        }
        return Fe;
    }
}