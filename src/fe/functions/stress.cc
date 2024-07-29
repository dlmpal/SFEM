#include "stress.h"
#include "../elements/solid/linear_elasticity.h"
#include "../elements/thermal/heat_conduction.h"

namespace sfem::fe::function
{
    //=============================================================================
    Stress2D::Stress2D()
        : Function(3)
    {
    }
    //=============================================================================
    la::DenseMatrix Stress2D::operator()(const FiniteElement &elem,
                                         const FEData &data,
                                         const std::vector<Scalar> &xpts,
                                         const std::vector<Scalar> &u) const
    {
        la::DenseMatrix stress(3, 1);

        if (elem.name() == "LinearElasticity2D")
        {
            auto solid_elem = dynamic_cast<const solid::LinearElasticity2D *>(&elem);
            auto constitutive = solid_elem->constitutive();
            stress = constitutive.eval_stress(data.dNdX, u);
        }

        return stress;
    }
    //=============================================================================
    VonMisesStress2D::VonMisesStress2D()
        : Function(1)
    {
    }
    //=============================================================================
    la::DenseMatrix VonMisesStress2D::operator()(const FiniteElement &elem,
                                                 const FEData &data,
                                                 const std::vector<Scalar> &xpts,
                                                 const std::vector<Scalar> &u) const
    {
        la::DenseMatrix stress_vm(1, 1);

        if (elem.name() == "LinearElasticity2D")
        {
            auto solid_elem = dynamic_cast<const solid::LinearElasticity2D *>(&elem);
            auto constitutive = solid_elem->constitutive();
            stress_vm.insert(0, 0, constitutive.eval_vm_stress(data.dNdX, u));
        }

        return stress_vm;
    }
    //=============================================================================
    Stress3D::Stress3D()
        : Function(6)
    {
    }
    //=============================================================================
    la::DenseMatrix Stress3D::operator()(const FiniteElement &elem,
                                         const FEData &data,
                                         const std::vector<Scalar> &xpts,
                                         const std::vector<Scalar> &u) const
    {
        la::DenseMatrix stress(6, 1);

        if (elem.name() == "LinearElasticity3D")
        {
            auto solid_elem = dynamic_cast<const solid::LinearElasticity3D *>(&elem);
            auto constitutive = solid_elem->constitutive();
            stress = constitutive.eval_stress(data.dNdX, u);
        }

        return stress;
    }
    //=============================================================================
    VonMisesStress3D::VonMisesStress3D()
        : Function(1)
    {
    }
    //=============================================================================
    la::DenseMatrix VonMisesStress3D::operator()(const FiniteElement &elem,
                                                 const FEData &data,
                                                 const std::vector<Scalar> &xpts,
                                                 const std::vector<Scalar> &u) const
    {
        la::DenseMatrix stress_vm(1, 1);

        if (elem.name() == "LinearElasticity3D")
        {
            auto solid_elem = dynamic_cast<const solid::LinearElasticity3D *>(&elem);
            auto constitutive = solid_elem->constitutive();
            stress_vm.insert(0, 0, constitutive.eval_vm_stress(data.dNdX, u));
        }

        return stress_vm;
    }
}