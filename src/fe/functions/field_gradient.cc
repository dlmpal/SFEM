#include "field_gradient.h"
#include "../finite_element.h"

namespace sfem::fe::function
{
    //=============================================================================
    FieldGradient::FieldGradient(const mesh::Field &field)
        : Function(field.n_vars() * 3)
    {
    }
    //=============================================================================
    la::DenseMatrix FieldGradient::operator()(const fe::FiniteElement &elem,
                                              const FEData &data,
                                              const std::vector<Scalar> &xpts,
                                              const std::vector<Scalar> &u,
                                              Scalar time) const
    {
        int n_vars = elem.n_vars();
        int n_nodes = elem.basis()->n_nodes();
        la::DenseMatrix grad(n_vars * 3, 1);

        for (int i = 0; i < n_vars; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < n_nodes; k++)
                {
                    grad.add(i * 3 + j, 0, data.dNdX[k * 3 + j] * u[k * n_vars + i]);
                }
            }
        }

        return grad;
    }
}