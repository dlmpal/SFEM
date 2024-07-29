#include "pressure_load.h"

namespace sfem::fe::solid
{
    //=============================================================================
    PressureLoad2D::PressureLoad2D(mesh::Cell cell, mesh::Field &P, Scalar thick)
        : FiniteElement("PressureLoad2D", 2, 2, cell), P_(P), thick_(thick)
    {
    }
    //=============================================================================
    la::DenseMatrix PressureLoad2D::evaluate_load_vector(const FEData &data,
                                                         const std::vector<Scalar> &xpts,
                                                         const std::vector<Scalar> &u) const
    {
        auto p = P_.get_cell_values(cell_);
        auto normal = cell_.face_normal(-1, xpts);
        la::DenseMatrix Fe(n_dof(), 1);
        for (int i = 0; i < basis_->n_nodes(); i++)
        {
            for (int j = 0; j < basis_->n_nodes(); j++)
            {
                Fe.add(i * 2 + 0, 0, -thick_ * p[j] * data.N[j] * normal.x[0] * data.N[i]);
                Fe.add(i * 2 + 1, 0, -thick_ * p[j] * data.N[j] * normal.x[1] * data.N[i]);
            }
        }
        return Fe;
    }
    //=============================================================================
    PressureLoad3D::PressureLoad3D(mesh::Cell cell, mesh::Field &P)
        : FiniteElement("PressureLoad3D", 3, 3, cell), P_(P)
    {
    }
    //=============================================================================
    la::DenseMatrix PressureLoad3D::evaluate_load_vector(const FEData &data,
                                                         const std::vector<Scalar> &xpts,
                                                         const std::vector<Scalar> &u) const
    {
        auto p = P_.get_cell_values(cell_);
        auto normal = cell_.face_normal(-1, xpts);
        la::DenseMatrix Fe(n_dof(), 1);
        for (int i = 0; i < basis_->n_nodes(); i++)
        {
            for (int j = 0; j < basis_->n_nodes(); j++)
            {
                Fe.add(i * 3 + 0, 0, -p[j] * data.N[j] * normal.x[0] * data.N[i]);
                Fe.add(i * 3 + 1, 0, -p[j] * data.N[j] * normal.x[1] * data.N[i]);
                Fe.add(i * 3 + 2, 0, -p[j] * data.N[j] * normal.x[2] * data.N[i]);
            }
        }
        return Fe;
    }
}