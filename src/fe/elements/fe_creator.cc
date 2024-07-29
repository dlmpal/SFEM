#include "fe_creator.h"

namespace sfem::fe
{
    std::shared_ptr<FiniteElement> FECreator::CreateSolidElement(const std::string &name,
                                                                 const mesh::Cell &cell,
                                                                 ThermoElasticConstitutive &constitutive) const
    {
        std::shared_ptr<FiniteElement> elem;

        if (name == "LinearElasticity2D")
        {
            elem = std::make_shared<solid::LinearElasticity2D>(cell, dynamic_cast<ThermoElasticPlaneConstitutive &>(constitutive));
        }
        else if (name == "LinearElasticity3D")
        {
            elem = std::make_shared<solid::LinearElasticity3D>(cell, dynamic_cast<ThermoElasticSolidConstitutive &>(constitutive));
        }
        else if (name == "HeatConduction2D")
        {
            elem = std::make_shared<thermal::HeatConduction2D>(cell, dynamic_cast<ThermoElasticPlaneConstitutive &>(constitutive));
        }
        else if (name == "HeatConduction3D")
        {
            elem = std::make_shared<thermal::HeatConduction3D>(cell, dynamic_cast<ThermoElasticSolidConstitutive &>(constitutive));
        }
        else
        {
            // error
        }

        return elem;
    }
}