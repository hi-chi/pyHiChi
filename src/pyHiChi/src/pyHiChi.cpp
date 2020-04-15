// Python-Interface.cpp : Defines exported functions for the dll-file.
//

#include <algorithm>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <pybind11/operators.h>

#include "pyGrid.h"

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


#define SET_METHODS_FOR_PY_GRID_FOR_FIELD_CONFIGURATIONS(pyGridType)     \
    .def("set", &pyGridType::setFieldConfiguration<NullField>)           \
    .def("set", &pyGridType::setFieldConfiguration<TightFocusingField>)    


#define SET_METHODS_FOR_PY_GRID(pyGridType)                              \
    .def(py::init<FP3, FP, FP3, FP3>())                                  \
    .def("getFields", &pyGridType::getFields)                            \
    .def("getJ", &pyGridType::getJ)                                      \
    .def("getE", &pyGridType::getE)                                      \
    .def("getB", &pyGridType::getB)                                      \
    .def("setJ", &pyGridType::setJ)                                      \
    .def("setE", &pyGridType::setE)                                      \
    .def("setB", &pyGridType::setB)                                      \
    .def("setJ", &pyGridType::pySetJ)                                    \
    .def("setE", &pyGridType::pySetE)                                    \
    .def("setB", &pyGridType::pySetB)                                    \
    .def("setJ", &pyGridType::setJxyz)                                   \
    .def("setE", &pyGridType::setExyz)                                   \
    .def("setB", &pyGridType::setBxyz)                                   \
    .def("setJ", &pyGridType::pySetJxyz)                                 \
    .def("setE", &pyGridType::pySetExyz)                                 \
    .def("setB", &pyGridType::pySetBxyz)                                 \
    .def("setJt", &pyGridType::setJxyzt)                                 \
    .def("setEt", &pyGridType::setExyzt)                                 \
    .def("setBt", &pyGridType::setBxyzt)                                 \
    SET_METHODS_FOR_PY_GRID_FOR_FIELD_CONFIGURATIONS(pyGridType)         \
    .def("analytical", &pyGridType::setAnalytical)                       \
    .def("setTime", &pyGridType::setTime)                                \
    .def("getTime", &pyGridType::getTime)


namespace py = pybind11;
using namespace pfc;


std::vector<ParticleType> ParticleInfo::typesVector = { {constants::electronMass, constants::electronCharge},//electron
                                    {constants::electronMass, -constants::electronCharge},//positron
                                    {constants::protonMass, 0.0},//proton
                                    {constants::electronMass, 0.0 } };//photon
const ParticleType* ParticleInfo::types = &ParticleInfo::typesVector[0];
short ParticleInfo::numTypes = sizeParticleTypes;

PYBIND11_MODULE(pyHiChi, object) {

    object.attr("pi") = constants::pi;
    object.attr("c") = constants::c;
    object.attr("lightVelocity") = constants::lightVelocity;
    object.attr("electronCharge") = constants::electronCharge;
    object.attr("electronMass") = constants::electronMass;
    object.attr("protonMass") = constants::protonMass;
    object.attr("planck") = constants::planck;
    object.attr("eV") = constants::eV;
    object.attr("meV") = constants::meV;

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

    py::class_<ValueField>(object, "field")
        .def(py::init<FP3, FP3>())
        .def(py::init<FP, FP, FP, FP, FP, FP>())
        .def("getE", &ValueField::getE)
        .def("setE", &ValueField::setE)
        .def("getB", &ValueField::getB)
        .def("setB", &ValueField::setB)
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
            return &(arr[i]);
        })
        .def("__setitem__", [](Ensemble3d &arr, size_t i, ParticleArray3d v) {
            if (i >= sizeParticleTypes) throw py::index_error();
            arr[i] = v;
        })
        .def("__getitem__", [](Ensemble3d& arr, string& name) {
            if (std::find(particleNames.begin(), particleNames.end(), name) == particleNames.end())
                throw py::index_error();
            return &(arr[name]);
        })
        .def("__setitem__", [](Ensemble3d &arr, string& name, ParticleArray3d v) {
            if (std::find(particleNames.begin(), particleNames.end(), name) == particleNames.end())
                throw py::index_error();
            arr[name] = v;
        })
        ;


    py::class_<BorisPusher>(object, "BorisPusher")
        .def(py::init<>())
        .def("__call__", (void (BorisPusher::*)(ParticleProxy3d*, ValueField, FP)) &BorisPusher::operator())
        .def("__call__", (void (BorisPusher::*)(Particle3d*, ValueField, FP)) &BorisPusher::operator())
        .def("__call__", (void (BorisPusher::*)(ParticleArray3d*, std::vector<ValueField>, FP)) &BorisPusher::operator())
        ;

    py::class_<RadiationReaction>(object, "RadiationReaction")
        .def(py::init<>())
        .def("__call__", (void (RadiationReaction::*)(ParticleProxy3d*, ValueField, FP)) &RadiationReaction::operator())
        .def("__call__", (void (RadiationReaction::*)(Particle3d*, ValueField, FP)) &RadiationReaction::operator())
        .def("__call__", (void (RadiationReaction::*)(ParticleArray3d*, std::vector<ValueField>, FP)) &RadiationReaction::operator())
        ;

    object.def("simpleThinning", &Thinning<ParticleArray3d>::simple);
    object.def("levelingThinning", &Thinning<ParticleArray3d>::leveling);
    object.def("numberConservativeThinning", &Thinning<ParticleArray3d>::numberConservative);
    object.def("energyConservativeThinning", &Thinning<ParticleArray3d>::energyConservative);
    object.def("kMeansMergining", &Merging<ParticleArray3d>::merge_with_kmeans);

    py::class_<pyYeeGrid>(object, "YeeGrid")
        SET_METHODS_FOR_PY_GRID(pyYeeGrid)
        ;

    py::class_<FDTD>(object, "FDTD")
        .def(py::init<pyYeeGrid*>())
        .def(py::init<pyYeeGridMapping*>())
        .def("setPML", &FDTD::setPML)
        .def("setBC", &FDTD::setFieldGenerator)
        .def("updateHalfB", &FDTD::updateHalfB)
        .def("updateE", &FDTD::updateE)
        .def("updateFields", &FDTD::updateFields)
        .def("setTimeStep", &FDTD::setTimeStep)
        ;

    py::class_<PeriodicalFieldGeneratorYee>(object, "PeriodicalBC")
        .def(py::init<RealFieldSolver<YeeGridType>*>())
        .def(py::init<FDTD*>())
        ;

    py::class_<pyPSTDGrid>(object, "PSTDGrid")
        SET_METHODS_FOR_PY_GRID(pyPSTDGrid)
        ;

    py::class_<PSTD>(object, "PSTD")
        .def(py::init<pyPSTDGrid*>())
        .def(py::init<pyPSTDGridMapping*>())
        .def("setPML", &PSTD::setPML)
        .def("updateFields", &PSTD::updateFields)
        .def("setTimeStep", &PSTD::setTimeStep)
        ;

    py::class_<pyPSATDGrid>(object, "PSATDGrid")
        SET_METHODS_FOR_PY_GRID(pyPSATDGrid)
        ;

    py::class_<PSATD>(object, "PSATD")
        .def(py::init<pyPSATDGrid*>())
        .def(py::init<pyPSATDGridMapping*>())
        .def("setPML", &PSATD::setPML)
        .def("updateFields", &PSATD::updateFields)
        .def("setTimeStep", &PSATD::setTimeStep)
        .def("convertFieldsPoissonEquation", &PSATD::convertFieldsPoissonEquation)
        ;
    
    py::class_<PSATDPoisson>(object, "PSATDPoisson")
        .def(py::init<pyPSATDGrid*>())
        .def(py::init<pyPSATDGridMapping*>())
        .def("setPML", &PSATDPoisson::setPML)
        .def("updateFields", &PSATDPoisson::updateFields)
        .def("setTimeStep", &PSATDPoisson::setTimeStep)
        .def("convertFieldsPoissonEquation", &PSATDPoisson::convertFieldsPoissonEquation)
        ;

    py::class_<pyPSATDTimeStraggeredGrid>(object, "PSATDTimeStraggeredGrid")
        SET_METHODS_FOR_PY_GRID(pyPSATDTimeStraggeredGrid)
        ;

    py::class_<PSATDTimeStraggered>(object, "PSATDTimeStraggered")
        .def(py::init<pyPSATDTimeStraggeredGrid*>())
        .def(py::init<pyPSATDTimeStraggeredGridMapping*>())
        .def("setPML", &PSATDTimeStraggered::setPML)
        .def("updateFields", &PSATDTimeStraggered::updateFields)
        .def("setTimeStep", &PSATDTimeStraggered::setTimeStep)
        .def("convertFieldsPoissonEquation", &PSATDTimeStraggered::convertFieldsPoissonEquation)
        ;
    
    py::class_<PSATDTimeStraggeredPoisson>(object, "PSATDTimeStraggeredPoisson")
        .def(py::init<pyPSATDTimeStraggeredGrid*>())
        .def(py::init<pyPSATDTimeStraggeredGridMapping*>())
        .def("setPML", &PSATDTimeStraggeredPoisson::setPML)
        .def("updateFields", &PSATDTimeStraggeredPoisson::updateFields)
        .def("setTimeStep", &PSATDTimeStraggeredPoisson::setTimeStep)
        .def("convertFieldsPoissonEquation", &PSATDTimeStraggeredPoisson::convertFieldsPoissonEquation)
        ;

    py::class_<ScalarQED_AEG_only_electron>(object, "QED")
        .def(py::init<>())
        .def("processParticles", &ScalarQED_AEG_only_electron::processParticles)
        .def("processParticlesNIter", &ScalarQED_AEG_only_electron::processParticlesNIter)
        ;

    py::class_<Mapping>(object, "Mapping")
        ;

    py::class_<pyYeeGridMapping>(object, "YeeGridMapping")
        .def(py::init<pyYeeGrid*>())
        SET_METHODS_FOR_PY_GRID(pyYeeGridMapping)
        .def("setMapping", &pyYeeGridMapping::setMapping)
        ;

    py::class_<pyPSTDGridMapping>(object, "PSTDGridMapping")
        .def(py::init<pyPSTDGrid*>())
        SET_METHODS_FOR_PY_GRID(pyPSTDGridMapping)
        .def("setMapping", &pyPSTDGridMapping::setMapping)
        ;

    py::class_<pyPSATDGridMapping>(object, "PSATDGridMapping")
        .def(py::init<pyPSATDGrid*>())
        SET_METHODS_FOR_PY_GRID(pyPSATDGridMapping)
        .def("setMapping", &pyPSATDGridMapping::setMapping)
        ;

    py::class_<pyPSATDTimeStraggeredGridMapping>(object, "PSATDTimeStraggeredGridMapping")
        .def(py::init<pyPSATDTimeStraggeredGrid*>())
        SET_METHODS_FOR_PY_GRID(pyPSATDTimeStraggeredGridMapping)
        .def("setMapping", &pyPSATDTimeStraggeredGridMapping::setMapping)
        ;

    // mappings

    py::class_<IdentityMapping, Mapping>(object, "IdentityMapping")
        .def(py::init<const FP3&, const FP3&>())
        .def("getDirectCoords", &IdentityMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &IdentityMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<PeriodicalMapping, Mapping>(object, "PeriodicalMapping")
        .def(py::init<Coordinate, FP, FP>())
        .def("getDirectCoords", &PeriodicalMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &PeriodicalMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<RotationMapping, Mapping>(object, "RotationMapping")
        .def(py::init<Coordinate, FP>())
        .def("getDirectCoords", &RotationMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &RotationMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<ScaleMapping, Mapping>(object, "ScaleMapping")
        .def(py::init<Coordinate, FP>())
        .def("getDirectCoords", &ScaleMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &ScaleMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<ShiftMapping, Mapping>(object, "ShiftMapping")
        .def(py::init<FP3>())
        .def("getDirectCoords", &ShiftMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &ShiftMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<TightFocusingMapping, Mapping>(object, "TightFocusingMapping")
        .def(py::init<FP, FP, FP>())
        .def("getDirectCoords", &TightFocusingMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getInverseCoords", &TightFocusingMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("getxMin", &TightFocusingMapping::getxMin)
        .def("getxMax", &TightFocusingMapping::getxMax)
        .def("setIfCut", &TightFocusingMapping::setIfCut)
        ;


    // field configurations

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
