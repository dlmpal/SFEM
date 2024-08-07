#include "inertial_load.h"

namespace sfem::fe::solid
{
    //=============================================================================
    InertialLoad2D::InertialLoad2D(const mesh::Cell cell,
                                   constitutive::ThermoElasticPlaneConstitutive &constitutive,
                                   const std::array<Scalar, 2> &g)
        : FiniteElement("InertialLoad2D", 2, 2, cell),
          constitutive_(constitutive),
          g_(g)
    {
    }
    //=============================================================================
    la::DenseMatrix InertialLoad2D::evaluate_load_vector(const FEData &data,
                                                         const std::vector<Scalar> &xpts,
                                                         const std::vector<Scalar> &u,
                                                         Scalar time) const
    {
        Scalar thick = constitutive_.thick();
        Scalar rho = constitutive_.prop().rho;
        la::DenseMatrix Fe(n_dof(), 1);
        for (int i = 0; i < basis_->n_nodes(); i++)
        {
            Fe.add(i * 2 + 0, 0, thick * rho * g_[0] * data.N[i]);
            Fe.add(i * 2 + 1, 0, thick * rho * g_[1] * data.N[i]);
        }
        return Fe;
    }
    //=============================================================================
    InertialLoad3D::InertialLoad3D(const mesh::Cell cell,
                                   constitutive::ThermoElasticSolidConstitutive &constitutive,
                                   const std::array<Scalar, 3> &g)
        : FiniteElement("InertialLoad3D", 3, 3, cell),
          constitutive_(constitutive),
          g_(g)
    {
    }
    //=============================================================================
    la::DenseMatrix InertialLoad3D::evaluate_load_vector(const FEData &data,
                                                         const std::vector<Scalar> &xpts,
                                                         const std::vector<Scalar> &u,
                                                         Scalar time) const
    {
        Scalar rho = constitutive_.prop().rho;
        la::DenseMatrix Fe(n_dof(), 1);
        for (int i = 0; i < basis_->n_nodes(); i++)
        {
            Fe.add(i * 3 + 0, 0, rho * g_[0] * data.N[i]);
            Fe.add(i * 3 + 1, 0, rho * g_[1] * data.N[i]);
            Fe.add(i * 3 + 2, 0, rho * g_[2] * data.N[i]);
        }
        return Fe;
    }
}