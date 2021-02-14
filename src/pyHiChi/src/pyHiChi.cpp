// Python-Interface.cpp : Defines exported functions for the dll-file.
//

#include <algorithm>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <pybind11/operators.h>

#include "pyField.h"

#include "Constants.h"
#include "Dimension.h"
#include "Ensemble.h"
#include "Fdtd.h"
#include "FieldGenerator.h"
#include "FieldValue.h"
#include "Merging.h"
#include "Particle.h"
#include "ParticleArray.h"
#include "ParticleTypes.h"
#include "Pstd.h"
#include "Psatd.h"
#include "Pusher.h"
#include "QED_AEG.h"
#include "Vectors.h"
#include "Thinning.h"
#include "Enums.h"
#include "Mapping.h"
#include "FieldConfiguration.h"


#define SET_FIELD_CONFIGURATIONS_GRID_METHODS(pyFieldType)                \
    .def("set", &pyFieldType::setFieldConfiguration<NullField>)           \
    .def("set", &pyFieldType::setFieldConfiguration<TightFocusingField>)    


#define SET_GRID_METHODS(pyFieldType)                                     \
    .def("set_J", &pyFieldType::setJ)                                      \
    .def("set_E", &pyFieldType::setE)                                      \
    .def("set_B", &pyFieldType::setB)                                      \
    .def("set_J", &pyFieldType::pySetJ)                                    \
    .def("set_E", &pyFieldType::pySetE)                                    \
    .def("set_B", &pyFieldType::pySetB)                                    \
    .def("set_J", &pyFieldType::setJxyz)                                   \
    .def("set_E", &pyFieldType::setExyz)                                   \
    .def("set_B", &pyFieldType::setBxyz)                                   \
    .def("set_J", &pyFieldType::pySetJxyz)                                 \
    .def("set_E", &pyFieldType::pySetExyz)                                 \
    .def("set_B", &pyFieldType::pySetBxyz)                                 \
    .def("set_J", &pyFieldType::setJxyzt)                                 \
    .def("set_E", &pyFieldType::setExyzt)                                 \
    .def("set_B", &pyFieldType::setBxyzt)                                 \
    SET_FIELD_CONFIGURATIONS_GRID_METHODS(pyFieldType)                    \
    .def("analytical", &pyFieldType::setAnalytical)                       \
    .def("set_time", &pyFieldType::setTime)                                \
    .def("get_time", &pyFieldType::getTime)                                \
    .def("refresh", &pyFieldType::refresh)


#define SET_FIELD_SOLVER_METHODS(pyFieldType)                             \
    .def("set_PML", &pyFieldType::setPML)                                  \
    .def("set_BC", &pyFieldType::setFieldGenerator)                        \
    .def("convert_fields_poisson_equation", &pyFieldType::convertFieldsPoissonEquation)  \
    .def("change_time_step", &pyFieldType::changeTimeStep)


#define SET_PY_FIELD_BASE_METHODS(pyFieldType)                            \
    .def("apply_mapping", [](std::shared_ptr<pyFieldType> self,            \
        std::shared_ptr<Mapping> mapping) {                               \
        return self->applyMapping(                                        \
            std::static_pointer_cast<pyFieldBase>(self), mapping          \
            );                                                            \
    })                                                                    \
    .def("__add__", [](std::shared_ptr<pyFieldType> self,                 \
        std::shared_ptr<pyFieldBase> other) {                             \
        return std::make_shared<pySumField>(                              \
            std::static_pointer_cast<pyFieldBase>(self), other            \
            );                                                            \
    }, py::is_operator())                                                 \
    .def("__mul__", [](std::shared_ptr<pyFieldType> self, FP factor) {    \
        return std::make_shared<pyMulField>(                              \
            std::static_pointer_cast<pyFieldBase>(self), factor           \
            );                                                            \
    }, py::is_operator())                                                 \
    .def("__rmul__", [](std::shared_ptr<pyFieldType> self, FP factor) {   \
        return std::make_shared<pyMulField>(                              \
            std::static_pointer_cast<pyFieldBase>(self), factor           \
            );                                                            \
    }, py::is_operator())


namespace py = pybind11;
using namespace pfc;


std::vector<ParticleType> ParticleInfo::typesVector = { {constants::electronMass, constants::electronCharge},//electron
                                    {constants::electronMass, -constants::electronCharge},//positron
                                    {constants::protonMass, 0.0},//proton
                                    {constants::electronMass, 0.0 } };//photon
const ParticleType* ParticleInfo::types = &ParticleInfo::typesVector[0];
short ParticleInfo::numTypes = sizeParticleTypes;

PYBIND11_MODULE(pyHiChi, object) {

    // ------------------- constants -------------------

    object.attr("pi") = constants::pi;
    object.attr("c") = constants::c;
    object.attr("LIGHT_VELOCITY") = constants::lightVelocity;
    object.attr("ELECTRON_CHARGE") = constants::electronCharge;
    object.attr("ELECTRON_MASS") = constants::electronMass;
    object.attr("PROTON_MASS") = constants::protonMass;
    object.attr("PLANCK") = constants::planck;
    object.attr("eV") = constants::eV;
    object.attr("meV") = constants::meV;

    // ------------------- auxulary structures -------------------

    py::enum_<Coordinate>(object, "Axis")
        .value("X", Coordinate::x)
        .value("Y", Coordinate::y)
        .value("Z", Coordinate::z)
        .export_values()
        ;

    py::class_<FP3>(object, "Vector3d")
        .def(py::init<>())
        .def(py::init<FP, FP, FP>())
        .def("volume", &FP3::volume)
        .def("norm", &FP3::norm)
        .def("norm2", &FP3::norm2)
        .def("normalize", &FP3::normalize)
        .def("__str__", &FP3::toString)

        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(-py::self)
        .def(FP() * py::self)
        .def(py::self * FP())
        .def(py::self *= FP())
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def(py::self / py::self)
        .def(py::self /= py::self)
        .def(py::self / FP())
        .def(py::self /= FP())

        .def_readwrite("x", &FP3::x)
        .def_readwrite("y", &FP3::y)
        .def_readwrite("z", &FP3::z)
        ;

    object.def("cross", (const Vector3<FP> (*)(const Vector3Proxy<FP>&, const Vector3Proxy<FP>&)) cross);
    object.def("cross", (FP3(*)(const FP3&, const FP3&)) cross);
    object.def("dot", (FP(*)(const Vector3Proxy<FP>&, const Vector3Proxy<FP>&)) dot);
    object.def("dot", (FP(*)(const FP3&, const FP3&)) dot);

    py::class_<ValueField>(object, "Field")
        .def(py::init<FP3, FP3>())
        .def(py::init<FP, FP, FP, FP, FP, FP>())
        .def("get_E", &ValueField::getE)
        .def("set_E", &ValueField::setE)
        .def("get_B", &ValueField::getB)
        .def("set_B", &ValueField::setB)
        ;

    // ------------------- particles -------------------

    py::class_<ParticleProxy3d>(object, "ParticleProxy")
        .def(py::init<Particle3d&>())
        .def(py::init<ParticleProxy3d&>())
        .def("get_position", &ParticleProxy3d::getPosition)
        .def("set_position", &ParticleProxy3d::setPosition)
        .def("get_momentum", &ParticleProxy3d::getMomentum)
        .def("set_momentum", &ParticleProxy3d::setMomentum)
        .def("get_velocity", &ParticleProxy3d::getVelocity)
        .def("set_velocity", &ParticleProxy3d::setVelocity)
        .def("get_weight", &ParticleProxy3d::getWeight)
        .def("set_weight", &ParticleProxy3d::setWeight)
        .def("get_gamma", &ParticleProxy3d::getGamma)
        .def("get_mass", &ParticleProxy3d::getMass)
        .def("get_charge", &ParticleProxy3d::getCharge)
        .def("get_type", &ParticleProxy3d::getType)
        ;

    py::enum_<ParticleTypes>(object, "ParticleTypes")
        .value("ELECTRON", Electron)
        .value("POSITRON", Positron)
        .value("PROTON", Proton)
        .export_values();

    py::class_<Particle3d>(object, "Particle")
        .def(py::init<>())
        .def(py::init<FP3, FP3>())
        .def(py::init<FP3, FP3, FP, ParticleTypes>())
        .def("get_position", &Particle3d::getPosition)
        .def("set_position", &Particle3d::setPosition)
        .def("get_momentum", &Particle3d::getMomentum)
        .def("set_momentum", &Particle3d::setMomentum)
        .def("get_velocity", &Particle3d::getVelocity)
        .def("set_velocity", &Particle3d::setVelocity)
        .def("get_weight", &Particle3d::getWeight)
        .def("set_weight", &Particle3d::setWeight)
        .def("get_gamma", &Particle3d::getGamma)
        .def("get_mass", &Particle3d::getMass)
        .def("get_charge", &Particle3d::getCharge)
        .def("get_type", &Particle3d::getType)
        ;

    py::class_<ParticleArray3d>(object, "ParticleArray")
        .def(py::init<>())
        .def(py::init<ParticleTypes>())
        .def("add", &ParticleArray3d::pushBack)
        .def("get_type", &ParticleArray3d::getType)
        .def("size", &ParticleArray3d::size)
        .def("delete", (void (ParticleArray3d::*)(int)) &ParticleArray3d::deleteParticle)
        .def("delete", (void (ParticleArray3d::*)(ParticleArray3d::iterator&)) &ParticleArray3d::deleteParticle)
        .def("__getitem__", [](ParticleArray3d& arr, size_t i) {
        if (i >= arr.size()) throw py::index_error();
        return arr[i];
    })
        .def("__setitem__", [](ParticleArray3d &arr, size_t i, Particle3d v) {
        if (i >= arr.size()) throw py::index_error();
        arr[i] = v;
    })
        .def("__iter__", [](ParticleArray3d &pArray) { return py::make_iterator(pArray.begin(), pArray.end()); },
            py::keep_alive<0, 1>())
        ;

    py::class_<Ensemble3d>(object, "Ensemble")
        .def(py::init<>())
        .def(py::init<Ensemble3d>())
        .def("add", &Ensemble3d::addParticle)
        .def("size", &Ensemble3d::size)
        .def("__getitem__", [](Ensemble3d& arr, size_t i) {
        if (i >= sizeParticleTypes) throw py::index_error();
        return std::reference_wrapper<Ensemble3d::ParticleArray>(arr[i]);
    })
        .def("__setitem__", [](Ensemble3d &arr, size_t i, ParticleArray3d v) {
        if (i >= sizeParticleTypes) throw py::index_error();
        arr[i] = v;
    })
        .def("__getitem__", [](Ensemble3d& arr, string& name) {
        if (std::find(particleNames.begin(), particleNames.end(), name) == particleNames.end())
            throw py::index_error();
        return std::reference_wrapper<Ensemble3d::ParticleArray>(arr[name]);
    })
        .def("__setitem__", [](Ensemble3d &arr, string& name, ParticleArray3d v) {
        if (std::find(particleNames.begin(), particleNames.end(), name) == particleNames.end())
            throw py::index_error();
        arr[name] = v;
    })
        ;

    // ------------------- pushers -------------------

    py::class_<BorisPusher>(object, "BorisPusher")
        .def(py::init<>())
        .def("__call__", (void (BorisPusher::*)(ParticleProxy3d*, ValueField&, FP)) &BorisPusher::operator())
        .def("__call__", (void (BorisPusher::*)(Particle3d*, ValueField&, FP)) &BorisPusher::operator())
        .def("__call__", (void (BorisPusher::*)(ParticleArray3d*, std::vector<ValueField>&, FP)) &BorisPusher::operator())
        ;

    // ------------------- other particle modules -------------------

    py::class_<RadiationReaction>(object, "RadiationReaction")
        .def(py::init<>())
        .def("__call__", (void (RadiationReaction::*)(ParticleProxy3d*, ValueField&, FP)) &RadiationReaction::operator())
        .def("__call__", (void (RadiationReaction::*)(Particle3d*, ValueField&, FP)) &RadiationReaction::operator())
        .def("__call__", (void (RadiationReaction::*)(ParticleArray3d*, std::vector<ValueField>&, FP)) &RadiationReaction::operator())
        ;

    py::class_<ScalarQED_AEG_only_electron>(object, "QED")
        .def(py::init<>())
        .def("process_particles", &ScalarQED_AEG_only_electron::processParticles)
        .def("process_particles", &ScalarQED_AEG_only_electron::processParticlesNIter)
        ;

    // ------------------- thinnings -------------------

    object.def("simple_thinning", &Thinning<ParticleArray3d>::simple);
    object.def("leveling_thinning", &Thinning<ParticleArray3d>::leveling);
    object.def("number_conservative_thinning", &Thinning<ParticleArray3d>::numberConservative);
    object.def("energy_conservative_thinning", &Thinning<ParticleArray3d>::energyConservative);
    object.def("k_means_mergining", &Merging<ParticleArray3d>::merge_with_kmeans);

    // ------------------- mappings -------------------

    py::class_<Mapping, std::shared_ptr<Mapping>> pyMapping(object, "Mapping");

    py::class_<IdentityMapping, std::shared_ptr<IdentityMapping>>(object, "IdentityMapping", pyMapping)
        .def(py::init<const FP3&, const FP3&>())
        .def("get_direct_coords", &IdentityMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &IdentityMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<PeriodicalMapping, std::shared_ptr<PeriodicalMapping>>(object, "PeriodicalMapping", pyMapping)
        .def(py::init<Coordinate, FP, FP>())
        .def("get_direct_coords", &PeriodicalMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &PeriodicalMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;
    
    py::class_<RotationMapping, std::shared_ptr<RotationMapping>>(object, "RotationMapping", pyMapping)
        .def(py::init<Coordinate, FP>())
        .def("get_direct_coords", &RotationMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &RotationMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<ScaleMapping, std::shared_ptr<ScaleMapping>>(object, "ScaleMapping", pyMapping)
        .def(py::init<Coordinate, FP>())
        .def("get_direct_coords", &ScaleMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &ScaleMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<ShiftMapping, std::shared_ptr<ShiftMapping>>(object, "ShiftMapping", pyMapping)
        .def(py::init<FP3>())
        .def("get_direct_coords", &ShiftMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &ShiftMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<TightFocusingMapping, std::shared_ptr<TightFocusingMapping>>(object, "TightFocusingMapping", pyMapping)
        .def(py::init<FP, FP, FP>())
        .def(py::init<FP, FP, FP, Coordinate>())
        .def("get_direct_coords", &TightFocusingMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &TightFocusingMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_min_coord", &TightFocusingMapping::getMinCoord)
        .def("get_max_coord", &TightFocusingMapping::getMaxCoord)
        .def("if_perform_inverse_mapping", &TightFocusingMapping::setIfCut)
        ;

    // ------------------- py fields -------------------

    //sample analytical fields
    py::class_<pyGaussianBeamField>(object, "GaussianBeamField")
        .def(py::init<FP, FP>())
        .def("get_fields", &pyGaussianBeamField::getFields)
        ;
    
    // abstract class
    py::class_<pyFieldBase, std::shared_ptr<pyFieldBase>> pyClassFieldBase(object, "FieldBase");
    pyClassFieldBase.def("get_fields", &pyFieldBase::getFields)
        .def("get_J", &pyFieldBase::getJ)
        .def("get_E", &pyFieldBase::getE)
        .def("get_B", &pyFieldBase::getB)
        .def("update_fields", &pyFieldBase::updateFields)
        .def("advance", &pyFieldBase::advance)
        ;

    py::class_<pySumField, std::shared_ptr<pySumField>>(
        object, "SumField", pyClassFieldBase)
        SET_PY_FIELD_BASE_METHODS(pySumField)
        ;

    py::class_<pyMulField, std::shared_ptr<pyMulField>>(
        object, "MulField", pyClassFieldBase)
        SET_PY_FIELD_BASE_METHODS(pyMulField)
        ;

    py::class_<pyYeeField, std::shared_ptr<pyYeeField>>(
        object, "YeeField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>())
        SET_GRID_METHODS(pyYeeField)
        SET_FIELD_SOLVER_METHODS(pyYeeField)
        SET_PY_FIELD_BASE_METHODS(pyYeeField)
        ;

    py::class_<pyPSTDField, std::shared_ptr<pyPSTDField>>(
        object, "PSTDField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>())
        SET_GRID_METHODS(pyPSTDField)
        SET_FIELD_SOLVER_METHODS(pyPSTDField)
        SET_PY_FIELD_BASE_METHODS(pyPSTDField)
        .def("set", &pyPSTDField::setEMField)
        .def("set", &pyPSTDField::pySetEMField)
        ;

    py::class_<pyPSATDField, std::shared_ptr<pyPSATDField>>(
        object, "PSATDField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>())
        SET_GRID_METHODS(pyPSATDField)
        SET_FIELD_SOLVER_METHODS(pyPSATDField)
        SET_PY_FIELD_BASE_METHODS(pyPSATDField)
        .def("set", &pyPSATDField::setEMField)
        .def("set", &pyPSATDField::pySetEMField)
        ;

    py::class_<pyPSATDPoissonField, std::shared_ptr<pyPSATDPoissonField>>(
        object, "PSATDPoissonField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>())
        SET_GRID_METHODS(pyPSATDPoissonField)
        SET_FIELD_SOLVER_METHODS(pyPSATDPoissonField)
        SET_PY_FIELD_BASE_METHODS(pyPSATDPoissonField)
        .def("set", &pyPSATDPoissonField::setEMField)
        .def("set", &pyPSATDPoissonField::pySetEMField)
        ;

    py::class_<pyPSATDTimeStraggeredField, std::shared_ptr<pyPSATDTimeStraggeredField>>(
        object, "PSATDTimeStraggeredField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>())
        SET_GRID_METHODS(pyPSATDTimeStraggeredField)
        SET_FIELD_SOLVER_METHODS(pyPSATDTimeStraggeredField)
        SET_PY_FIELD_BASE_METHODS(pyPSATDTimeStraggeredField)
        .def("set", &pyPSATDTimeStraggeredField::setEMField)
        .def("set", &pyPSATDTimeStraggeredField::pySetEMField)
        ;

    py::class_<pyPSATDTimeStraggeredPoissonField, std::shared_ptr<pyPSATDTimeStraggeredPoissonField>>(
        object, "PSATDTimeStraggeredPoissonField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>())
        SET_GRID_METHODS(pyPSATDTimeStraggeredPoissonField)
        SET_FIELD_SOLVER_METHODS(pyPSATDTimeStraggeredPoissonField)
        SET_PY_FIELD_BASE_METHODS(pyPSATDTimeStraggeredPoissonField)
        .def("set", &pyPSATDTimeStraggeredPoissonField::setEMField)
        .def("set", &pyPSATDTimeStraggeredPoissonField::pySetEMField)
        ;

    // ------------------- field generators -------------------

    py::class_<PeriodicalFieldGeneratorYee, std::shared_ptr<PeriodicalFieldGeneratorYee>>(object, "PeriodicalBC")
        //.def(py::init<RealFieldSolver<GridTypes::YeeGridType>*>())
        .def(py::init([](std::shared_ptr<pyYeeField> field) {
        return std::make_shared<PeriodicalFieldGeneratorYee>(field->getFieldSolver());
    }))
        ;

    // ------------------- field configurations -------------------

    py::class_<NullField>(object, "NullField")
        .def(py::init<>())
        .def("get_E", &NullField::getE)
        .def("get_B", &NullField::getB)
        ;

    py::class_<TightFocusingField>(object, "TightFocusingField")
        .def(py::init<FP, FP, FP, FP, FP, FP>())
        .def(py::init<FP, FP, FP, FP, FP, FP, FP3>())
        .def(py::init<FP, FP, FP, FP, FP, FP, FP3, FP>())
        .def("get_E", &TightFocusingField::getE)
        .def("get_B", &TightFocusingField::getB)
        ;

}
