#pragma once

#include "solid/sfem_solid_elements.h"
#include "thermal/sfem_thermal_elements.h"

namespace sfem::fe
{
    using namespace constitutive;

    class FECreator
    {
    public:
        FECreator() = default;

        std::shared_ptr<FiniteElement> CreateSolidElement(const std::string &name,
                                                          const mesh::Cell &cell,
                                                          ThermoElasticConstitutive &constitutive) const;
        // std ::shared_ptr<FiniteElement> CreatePressureLoadElement() const;
        // std::shared_ptr<FiniteElement> CreateThermalLoadElement() const;
        // std::shared_ptr<FiniteElement> CreateHeatFluxElement() const;
        // std::shared_ptr<FiniteElement> CreateHeatConvectionElement() const;
    };
}