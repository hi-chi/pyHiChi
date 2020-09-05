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
#include "Handler.h"
#include "Merging.h"
#include "Particle.h"
#include "ParticleArray.h"
#include "ParticleTypes.h"
#include "Pstd.h"
#include "Psatd.h"
#include "QED_AEG.h"
#include "Vectors.h"
#include "Thinning.h"
#include "Enums.h"
#include "Mapping.h"
#include "FieldConfiguration.h"


#define SET_METHODS_OF_PY_GRID_FIELD_CONFIGURATIONS(pyFieldType)          \
    .def("set", &pyFieldType::setFieldConfiguration<NullField>)           \
    .def("set", &pyFieldType::setFieldConfiguration<TightFocusingField>)    


#define SET_METHODS_OF_PY_GRID(pyFieldType)                               \
    .def("getFields", &pyFieldType::getFields)                            \
    .def("getJ", &pyFieldType::getJ)                                      \
    .def("getE", &pyFieldType::getE)                                      \
    .def("getB", &pyFieldType::getB)                                      \
    .def("setJ", &pyFieldType::setJ)                                      \
    .def("setE", &pyFieldType::setE)                                      \
    .def("setB", &pyFieldType::setB)                                      \
    .def("setJ", &pyFieldType::pySetJ)                                    \
    .def("setE", &pyFieldType::pySetE)                                    \
    .def("setB", &pyFieldType::pySetB)                                    \
    .def("setJ", &pyFieldType::setJxyz)                                   \
    .def("setE", &pyFieldType::setExyz)                                   \
    .def("setB", &pyFieldType::setBxyz)                                   \
    .def("setJ", &pyFieldType::pySetJxyz)                                 \
    .def("setE", &pyFieldType::pySetExyz)                                 \
    .def("setB", &pyFieldType::pySetBxyz)                                 \
    .def("setJt", &pyFieldType::setJxyzt)                                 \
    .def("setEt", &pyFieldType::setExyzt)                                 \
    .def("setBt", &pyFieldType::setBxyzt)                                 \
    SET_METHODS_OF_PY_GRID_FIELD_CONFIGURATIONS(pyFieldType)              \
    .def("analytical", &pyFieldType::setAnalytical)                       \
    .def("setTime", &pyFieldType::setTime)                                \
    .def("getTime", &pyFieldType::getTime)                                \
    .def("refresh", &pyFieldType::refresh)


#define SET_METHODS_OF_FIELD_SOLVER(pyFieldType)                          \
    .def("setPML", &pyFieldType::setPML)                                  \
    .def("setBC", &pyFieldType::setFieldGenerator)                        \
    .def("updateFields", &pyFieldType::updateFields)                      \
    .def("advance", &pyFieldType::advance)                                \
    .def("changeTimeStep", &pyFieldType::changeTimeStep)


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
    object.attr("lightVelocity") = constants::lightVelocity;
    object.attr("electronCharge") = constants::electronCharge;
    object.attr("electronMass") = constants::electronMass;
    object.attr("protonMass") = constants::protonMass;
    object.attr("planck") = constants::planck;
    object.attr("eV") = constants::eV;
    object.attr("meV") = constants::meV;

    // ------------------- auxulary structures -------------------

    py::enum_<Coordinate>(object, "Axis")
        .value("x", Coordinate::x)
        .value("y", Coordinate::y)
        .value("z", Coordinate::z)
        .export_values()
        ;

    py::class_<FP3>(object, "vector3d")
        .def(py::init<>())
        .def(py::init<FP, FP, FP>())
        .def("volume", &FP3::volume)
        .def("norm", &FP3::norm)
        .def("norm2", &FP3::norm2)
        .def("normalize", &FP3::normalize)
        .def("toString", &FP3::toString)
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
    object.def("cross", cross);
    object.def("dot", dot);

    py::class_<ValueField>(object, "field")
        .def(py::init<FP3, FP3>())
        .def(py::init<FP, FP, FP, FP, FP, FP>())
        .def("getE", &ValueField::getE)
        .def("setE", &ValueField::setE)
        .def("getB", &ValueField::getB)
        .def("setB", &ValueField::setB)
        ;

    // ------------------- particles -------------------

    py::class_<ParticleProxy3d>(object, "particleProxy")
        .def(py::init<Particle3d&>())
        .def(py::init<ParticleProxy3d&>())
        .def("getPosition", &ParticleProxy3d::getPosition)
        .def("setPosition", &ParticleProxy3d::setPosition)
        .def("getMomentum", &ParticleProxy3d::getMomentum)
        .def("setMomentum", &ParticleProxy3d::setMomentum)
        .def("getVelocity", &ParticleProxy3d::getVelocity)
        .def("setVelocity", &ParticleProxy3d::setVelocity)
        .def("getWeight", &ParticleProxy3d::getWeight)
        .def("setWeight", &ParticleProxy3d::setWeight)
        .def("getGamma", &ParticleProxy3d::getGamma)
        .def("getMass", &ParticleProxy3d::getMass)
        .def("getCharge", &ParticleProxy3d::getCharge)
        .def("getType", &ParticleProxy3d::getType)
        ;
    
    py::enum_<ParticleTypes>(object, "particleTypes")
        .value("Electron", Electron)
        .value("Positron", Positron)
        .value("Proton", Proton)
        .value("Photon", Photon)
        .export_values();

    py::class_<Particle3d>(object, "particle")
        .def(py::init<>())
        .def(py::init<FP3, FP3>())
        .def(py::init<FP3, FP3, FP, ParticleTypes>())
        .def("getPosition", &Particle3d::getPosition)
        .def("setPosition", &Particle3d::setPosition)
        .def("getMomentum", &Particle3d::getMomentum)
        .def("setMomentum", &Particle3d::setMomentum)
        .def("getVelocity", &Particle3d::getVelocity)
        .def("setVelocity", &Particle3d::setVelocity)
        .def("getWeight", &Particle3d::getWeight)
        .def("setWeight", &Particle3d::setWeight)
        .def("getGamma", &Particle3d::getGamma)
        .def("getMass", &Particle3d::getMass)
        .def("getCharge", &Particle3d::getCharge)
        .def("getType", &Particle3d::getType)
        ;
    
    py::class_<ParticleArray3d>(object, "particleArray")
        .def(py::init<>())
        .def(py::init<ParticleTypes>())
        .def("add", &ParticleArray3d::pushBack)
        .def("getType", &ParticleArray3d::getType)
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

    py::class_<Ensemble3d>(object, "ensemble")
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
        .def("__call__", (void (BorisPusher::*)(ParticleProxy3d*, ValueField, FP)) &BorisPusher::operator())
        .def("__call__", (void (BorisPusher::*)(Particle3d*, ValueField, FP)) &BorisPusher::operator())
        .def("__call__", (void (BorisPusher::*)(ParticleArray3d*, std::vector<ValueField>, FP)) &BorisPusher::operator())
        ;

    // ------------------- other particle modules -------------------

    py::class_<RadiationReaction>(object, "RadiationReaction")
        .def(py::init<>())
        .def("__call__", (void (RadiationReaction::*)(ParticleProxy3d*, ValueField, FP)) &RadiationReaction::operator())
        .def("__call__", (void (RadiationReaction::*)(Particle3d*, ValueField, FP)) &RadiationReaction::operator())
        .def("__call__", (void (RadiationReaction::*)(ParticleArray3d*, std::vector<ValueField>, FP)) &RadiationReaction::operator())
        ;

    py::class_<ScalarQED_AEG_only_electron>(object, "QED")
        .def(py::init<>())
        .def("processParticles", &ScalarQED_AEG_only_electron::processParticles)
        .def("processParticlesNIter", &ScalarQED_AEG_only_electron::processParticlesNIter)
        ;

    // ------------------- thinnings -------------------

    object.def("simpleThinning", &Thinning<ParticleArray3d>::simple);
    object.def("levelingThinning", &Thinning<ParticleArray3d>::leveling);
    object.def("numberConservativeThinning", &Thinning<ParticleArray3d>::numberConservative);
    object.def("energyConservativeThinning", &Thinning<ParticleArray3d>::energyConservative);
    object.def("kMeansMergining", &Merging<ParticleArray3d>::merge_with_kmeans);

    // ------------------- mappings -------------------

    py::class_<Mapping, std::shared_ptr<Mapping>> pyMapping(object, "Mapping");

    py::class_<IdentityMapping, std::shared_ptr<IdentityMapping>>(object, "IdentityMapping", pyMapping)
        .def(py::init<const FP3&, const FP3&>())
        .def("getDirectCoords", &IdentityMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &IdentityMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<PeriodicalMapping, std::shared_ptr<PeriodicalMapping>>(object, "PeriodicalMapping", pyMapping)
        .def(py::init<Coordinate, FP, FP>())
        .def("getDirectCoords", &PeriodicalMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &PeriodicalMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<RotationMapping, std::shared_ptr<RotationMapping>>(object, "RotationMapping", pyMapping)
        .def(py::init<Coordinate, FP>())
        .def("getDirectCoords", &RotationMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &RotationMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<ScaleMapping, std::shared_ptr<ScaleMapping>>(object, "ScaleMapping", pyMapping)
        .def(py::init<Coordinate, FP>())
        .def("getDirectCoords", &ScaleMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &ScaleMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<ShiftMapping, std::shared_ptr<ShiftMapping>>(object, "ShiftMapping", pyMapping)
        .def(py::init<FP3>())
        .def("getDirectCoords", &ShiftMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &ShiftMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<TightFocusingMapping, std::shared_ptr<TightFocusingMapping>>(object, "TightFocusingMapping", pyMapping)
        .def(py::init<FP, FP, FP>())
        .def(py::init<FP, FP, FP, Coordinate>())
        .def("getDirectCoords", &TightFocusingMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &TightFocusingMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getMinCoord", &TightFocusingMapping::getMinCoord)
        .def("getMaxCoord", &TightFocusingMapping::getMaxCoord)
        .def("setIfCut", &TightFocusingMapping::setIfCut)
        ;

    // ------------------- py fields -------------------

    py::class_<pyYeeField, std::shared_ptr<pyYeeField>>(object, "YeeField")
        .def(py::init<FP3, FP3, FP3, FP>())
        .def("applyMapping", [](std::shared_ptr<pyYeeField> self, std::shared_ptr<Mapping> mapping) {
        return pyYeeField::applyMapping(self, mapping);
    })
        SET_METHODS_OF_PY_GRID(pyYeeField)
        SET_METHODS_OF_FIELD_SOLVER(pyYeeField)
        ;

    py::class_<pyPSTDField, std::shared_ptr<pyPSTDField>>(object, "PSTDField")
        .def(py::init<FP3, FP3, FP3, FP>())
        .def("applyMapping", [](std::shared_ptr<pyPSTDField> self,
            std::shared_ptr<Mapping> mapping) {
        return pyPSTDField::applyMapping(self, mapping);
    })
        SET_METHODS_OF_PY_GRID(pyPSTDField)
        SET_METHODS_OF_FIELD_SOLVER(pyPSTDField)
        .def("set", &pyPSTDField::setEMField)
        .def("set", &pyPSTDField::pySetEMField)
        ;


    py::class_<pyPSATDField, std::shared_ptr<pyPSATDField>>(object, "PSATDField")
        .def(py::init<FP3, FP3, FP3, FP>())
        .def("applyMapping", [](std::shared_ptr<pyPSATDField> self,
            std::shared_ptr<Mapping> mapping) {
        return pyPSATDField::applyMapping(self, mapping);
    })
        SET_METHODS_OF_PY_GRID(pyPSATDField)
        SET_METHODS_OF_FIELD_SOLVER(pyPSATDField)
        .def("set", &pyPSATDField::setEMField)
        .def("set", &pyPSATDField::pySetEMField)
        ;


    py::class_<pyPSATDPoissonField, std::shared_ptr<pyPSATDPoissonField>>(object, "PSATDPoissonField")
        .def(py::init<FP3, FP3, FP3, FP>())
        .def("applyMapping", [](std::shared_ptr<pyPSATDPoissonField> self,
            std::shared_ptr<Mapping> mapping) {
        return pyPSATDPoissonField::applyMapping(self, mapping);
    })
        SET_METHODS_OF_PY_GRID(pyPSATDPoissonField)
        SET_METHODS_OF_FIELD_SOLVER(pyPSATDPoissonField)
        .def("set", &pyPSATDPoissonField::setEMField)
        .def("set", &pyPSATDPoissonField::pySetEMField)
        ;


    py::class_<pyPSATDTimeStraggeredField, std::shared_ptr<pyPSATDTimeStraggeredField>>(
        object, "PSATDTimeStraggeredField")
        .def(py::init<FP3, FP3, FP3, FP>())
        .def("applyMapping", [](std::shared_ptr<pyPSATDTimeStraggeredField> self,
            std::shared_ptr<Mapping> mapping) {
        return pyPSATDTimeStraggeredField::applyMapping(self, mapping);
    })
        SET_METHODS_OF_PY_GRID(pyPSATDTimeStraggeredField)
        SET_METHODS_OF_FIELD_SOLVER(pyPSATDTimeStraggeredField)
        .def("set", &pyPSATDTimeStraggeredField::setEMField)
        .def("set", &pyPSATDTimeStraggeredField::pySetEMField)
        ;


    py::class_<pyPSATDTimeStraggeredPoissonField, std::shared_ptr<pyPSATDTimeStraggeredPoissonField>>(
        object, "PSATDTimeStraggeredPoissonField")
        .def(py::init<FP3, FP3, FP3, FP>())
        .def("applyMapping", [](std::shared_ptr<pyPSATDTimeStraggeredPoissonField> self,
            std::shared_ptr<Mapping> mapping) {
        return pyPSATDTimeStraggeredPoissonField::applyMapping(self, mapping);
    })
        SET_METHODS_OF_PY_GRID(pyPSATDTimeStraggeredPoissonField)
        SET_METHODS_OF_FIELD_SOLVER(pyPSATDTimeStraggeredPoissonField)
        .def("set", &pyPSATDTimeStraggeredPoissonField::setEMField)
        .def("set", &pyPSATDTimeStraggeredPoissonField::pySetEMField)
        ;


    // ------------------- field generators -------------------

    py::class_<PeriodicalFieldGeneratorYee>(object, "PeriodicalBC")
        .def(py::init<RealFieldSolver<GridTypes::YeeGridType>*>())
        .def(py::init<FDTD*>())
        ;

    // ------------------- field configurations -------------------

    py::class_<NullField>(object, "NullField")
        .def(py::init<>())
        .def("getE", &NullField::getE)
        .def("getB", &NullField::getB)
        ;

    py::class_<TightFocusingField>(object, "TightFocusingField")
        .def(py::init<FP, FP, FP, FP, FP, FP, FP>())
        .def(py::init<FP, FP, FP, FP, FP, FP, FP, FP3>())
        .def(py::init<FP, FP, FP, FP, FP, FP, FP, FP3, FP>())
        .def("getE", &TightFocusingField::getE)
        .def("getB", &TightFocusingField::getB)
        ;

}
