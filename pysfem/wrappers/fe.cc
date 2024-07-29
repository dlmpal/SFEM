#include "sfem.h"
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

using namespace sfem;
using namespace fe;
namespace nb = nanobind;

namespace sfem_wrappers
{
    //=============================================================================
    void init_constitutive(nb::module_ &m);
    void init_thermal(nb::module_ &m);
    void init_solid(nb::module_ &m);
    void init_function(nb::module_ &m);
    void init_utils(nb::module_ &m);
    //=============================================================================
    void init_fe(nb::module_ &m)
    {
        // FEData
        nb::class_<FEData>(m, "FEData")
            .def_rw("qwt", &FEData::qwt)
            .def_rw("qpt", &FEData::qpt)
            .def_rw("detJ", &FEData::detJ)
            .def_rw("dXdxi", &FEData::dXdxi)
            .def_rw("dxidX", &FEData::dxidX)
            .def_rw("N", &FEData::N)
            .def_rw("dNdxi", &FEData::dNdxi)
            .def_rw("dNdX", &FEData::dNdX);

        // FEMatrixType
        nb::enum_<FEMatrixType>(m, "FEMatrixType")
            .value("mass", FEMatrixType::mass)
            .value("damping", FEMatrixType::damping)
            .value("stiffness", FEMatrixType::stiffness)
            .value("jacobian", FEMatrixType::jacobian);

        // FEVectorType
        nb::enum_<FEVectorType>(m, "FEVectorType")
            .value("load", FEVectorType::load)
            .value("residual", FEVectorType::residual);

        // FiniteElement
        nb::class_<FiniteElement>(m, "FiniteElement")
            .def(nb::init<const std::string &, int, int, mesh::Cell>())
            .def("name", &FiniteElement::name)
            .def("n_vars", &FiniteElement::n_vars)
            .def("physical_dim", &FiniteElement::physical_dim)
            .def("n_dof", &FiniteElement::n_dof)
            .def("cell", &FiniteElement::cell, nb::rv_policy::reference_internal)
            .def("basis", &FiniteElement::basis, nb::rv_policy::reference_internal)
            .def("transform_basis", &FiniteElement::transform_basis)
            .def("evaluate_mass_matrix", &FiniteElement::evaluate_mass_matrix)
            .def("evaluate_damping_matrix", &FiniteElement::evaluate_damping_matrix)
            .def("evaluate_stiff_matrix", &FiniteElement::evaluate_stiff_matrix)
            .def("evaluate_load_vector", &FiniteElement::evaluate_load_vector)
            .def("integrate_fe_matrix", &FiniteElement::integrate_fe_matrix)
            .def("integrate_fe_vector", &FiniteElement::integrate_fe_vector);

        // Constitutive
        nb::module_ constitutive = m.def_submodule("constitutive", "Constitutive laws");
        init_constitutive(constitutive);

        // Thermal elements
        // nb::module_ thermal = m.def_submodule("thermal", "Thermal elements");
        // init_thermal(thermal);

        // Solid elements
        nb::module_ solid = m.def_submodule("solid", "Solid elements");
        init_solid(solid);

        // FECreator
        nb::class_<fe::FECreator>(m, "FECreator")
            .def(nb::init<>())
            .def("create_solid_element", &FECreator::CreateSolidElement);

        // Functions
        nb::module_ function = m.def_submodule("function", "Finite element functions");
        init_function(function);

        // Utils
        init_utils(m);
    }
    //=============================================================================
    void init_constitutive(nb::module_ &m)
    {
        using namespace constitutive;

        // ThermoMechanicalProperties
        nb::class_<ThermoMechanicalProperties>(m, "ThermoMechanicalProprties")
            .def(nb::init<>())
            .def_rw("E", &ThermoMechanicalProperties::E, "Young's modulus")
            .def_rw("nu", &ThermoMechanicalProperties::nu, "Poisson's ratio")
            .def_rw("kappa", &ThermoMechanicalProperties::kappa, "Thermal conductivity")
            .def_rw("rho", &ThermoMechanicalProperties::rho, "Density")
            .def_rw("cp", &ThermoMechanicalProperties::cp, "Specific heat")
            .def_rw("alpha", &ThermoMechanicalProperties::alpha, "Heat expansion coefficient");

        // ThermoElasticConstitutive
        nb::class_<ThermoElasticConstitutive>(m, "ThermoElasticConstitutive")
            .def("prop", &ThermoElasticConstitutive::prop, nb::rv_policy::reference_internal);

        // ThermoElasticSolidConstitutive
        nb::class_<ThermoElasticSolidConstitutive, ThermoElasticConstitutive>(m, "ThermoElasticSolidConstitutive")
            .def(nb::init<ThermoMechanicalProperties &>());
    }
    //=============================================================================
    void init_thermal(nb::module_ &m)
    {
    }
    //=============================================================================
    void init_solid(nb::module_ &m)
    {
        using namespace solid;

        // LinearElasticity2D
        nb::class_<LinearElasticity2D, FiniteElement>(m, "LinearElasticity2D")
            .def(nb::init<mesh::Cell, ThermoElasticPlaneConstitutive &>())
            .def("constitutive", &LinearElasticity2D::constitutive)
            .def("constitutive", &LinearElasticity2D::constitutive)
            .def("add_inertial_load", &LinearElasticity2D::add_inertial_load)
            .def("add_thermal_load", &LinearElasticity2D::add_thermal_load);

        // LinearElasticity3D
        nb::class_<LinearElasticity3D, FiniteElement>(m, "LinearElasticity3D")
            .def("constitutive", &LinearElasticity3D::constitutive)
            .def(nb::init<mesh::Cell, ThermoElasticSolidConstitutive &>())
            .def("constitutive", &LinearElasticity3D::constitutive)
            .def("add_inertial_load", &LinearElasticity3D::add_inertial_load);
    }
    //=============================================================================
    void init_function(nb::module_ &m)
    {
        // Function
        nb::class_<Function>(m, "Function");

        using namespace function;

        // FieldGradient
        nb::class_<FieldGradient, Function>(m, "FieldGradient")
            .def(nb::init<const mesh::Field &>());

        // StructuralMass2D
        nb::class_<StructuralMass2D, Function>(m, "StructuralMass2D")
            .def(nb::init<>());

        // StructuralMass3D
        nb::class_<StructuralMass3D, Function>(m, "StructuralMass3D")
            .def(nb::init<>());

        // Stress2D
        nb::class_<Stress2D, Function>(m, "Stress2D")
            .def(nb::init<>());

        // Stress3D
        nb::class_<Stress3D, Function>(m, "Stress3D")
            .def(nb::init<>());

        // VonMisesStress2D
        nb::class_<VonMisesStress2D, Function>(m, "VonMisesStress2D")
            .def(nb::init<>());

        // VonMisesStress3D
        nb::class_<VonMisesStress3D, Function>(m, "VonMisesStress3D")
            .def(nb::init<>());
    }
    //=============================================================================
    void init_utils(nb::module_ &m)
    {
        // Assembly
        m.def("assemble_function", &assemble_function);
        m.def("assemble_vector", &assemble_vector);
        m.def("assemble_matrix", &assemble_matrix);

        // Project
        m.def("project_function", &project_function);
    }
}
