#include "thermal_load.h"

namespace sfem::fe::solid
{
    //=============================================================================
    ThermalLoad2D::ThermalLoad2D(mesh::Cell cell,
                                 constitutive::ThermoElasticPlaneConstitutive &constitutive,
                                 mesh::Field &T,
                                 Scalar T0)
        : FiniteElement("ThermalLoad2D", 2, 2, cell),
          constitutive_(constitutive),
          T_(T),
          T0_(T0)
    {
    }
    //=============================================================================
    la::DenseMatrix ThermalLoad2D::evaluate_load_vector(const FEData &data,
                                                        const std::vector<Scalar> &xpts,
                                                        const std::vector<Scalar> &u,
                                                        Scalar time) const
    {
        auto dT = T_.get_cell_values(cell_);
        for (int i = 0; i < basis_->n_nodes(); i++)
        {
            dT[i] -= T0_;
        }
        return constitutive_.eval_thermal_stress(data.N, data.dNdX, dT);
    }
}