#include "structural_mass.h"
#include "../elements/solid/linear_elasticity.h"
#include "../elements/thermal/heat_conduction.h"

namespace sfem::fe::function
{
    //=============================================================================
    StructuralMass2D::StructuralMass2D()
        : Function(1)
    {
    }
    //=============================================================================
    la::DenseMatrix StructuralMass2D::operator()(const FiniteElement &elem,
                                                 const FEData &data,
                                                 const std::vector<Scalar> &xpts,
                                                 const std::vector<Scalar> &u,
                                                 Scalar time) const
    {
        la::DenseMatrix mass(1, 1);

        if (elem.name() == "LinearElasticity2D")
        {
            auto solid_elem = dynamic_cast<const solid::LinearElasticity2D *>(&elem);
            auto constitutive = solid_elem->constitutive();
            mass.add(0, 0, constitutive.prop().rho * constitutive.thick());
        }
        else if (elem.name() == "HeatConduction2D")
        {
            auto thermal_elem = dynamic_cast<const thermal::HeatConduction2D *>(&elem);
            auto constitutive = thermal_elem->constitutive();
            mass.add(0, 0, constitutive.prop().rho * constitutive.thick());
        }

        return mass;
    }
    //=============================================================================
    StructuralMass3D::StructuralMass3D()
        : Function(1)
    {
    }
    //=============================================================================
    la::DenseMatrix StructuralMass3D::operator()(const FiniteElement &elem,
                                                 const FEData &data,
                                                 const std::vector<Scalar> &xpts,
                                                 const std::vector<Scalar> &u,
                                                 Scalar time) const
    {
        la::DenseMatrix mass(1, 1);

        if (elem.name() == "LinearElasticity3D")
        {
            auto solid_elem = dynamic_cast<const solid::LinearElasticity3D *>(&elem);
            auto constitutive = solid_elem->constitutive();
            mass.add(0, 0, constitutive.prop().rho);
        }
        else if (elem.name() == "HeatConduction3D")
        {
            auto thermal_elem = dynamic_cast<const thermal::HeatConduction3D *>(&elem);
            auto constitutive = thermal_elem->constitutive();
            mass.add(0, 0, constitutive.prop().rho);
        }

        return mass;
    }
}