// Python-Interface.cpp : Определяет экспортированные функции для приложения DLL.
//

#include <algorithm>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

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

	py::class_<FP3>(object, "vector3d")
		.def(py::init<>())
		.def(py::init<FP, FP, FP>())
		.def("volume", &FP3::volume)
		.def("norm", &FP3::norm)
		.def("norm2", &FP3::norm2)
        .def("normalize", &FP3::normalize)
		.def("toString", &FP3::toString)
		.def("__str__", &FP3::toString)

		.def_readwrite("x", &FP3::x)
		.def_readwrite("y", &FP3::y)
		.def_readwrite("z", &FP3::z)
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
		.def(py::init<FP3, FP, FP3, FP3>())
		.def("getFields", &pyYeeGrid::getFields)
		.def("getJ", &pyYeeGrid::getJ)
		.def("getE", &pyYeeGrid::getE)
		.def("getB", &pyYeeGrid::getB)
		.def("setJ", &pyYeeGrid::setJ)
		.def("setE", &pyYeeGrid::setE)
		.def("setB", &pyYeeGrid::setB)
		.def("setJ", &pyYeeGrid::pySetJ)
		.def("setE", &pyYeeGrid::pySetE)
		.def("setB", &pyYeeGrid::pySetB)
		.def("setJ", &pyYeeGrid::setJxyz)
		.def("setE", &pyYeeGrid::setExyz)
		.def("setB", &pyYeeGrid::setBxyz)
		.def("setJ", &pyYeeGrid::pySetJxyz)
		.def("setE", &pyYeeGrid::pySetExyz)
		.def("setB", &pyYeeGrid::pySetBxyz)
		.def("setJt", &pyYeeGrid::setJxyzt)
		.def("setEt", &pyYeeGrid::setExyzt)
		.def("setBt", &pyYeeGrid::setBxyzt)
		.def("analytical", &pyYeeGrid::setAnalytical)
		.def("setTime", &pyYeeGrid::setTime)
		;

	py::class_<FDTD>(object, "FDTD")
		.def(py::init<pyYeeGrid*>())
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
        .def(py::init<FP3, FP, FP3, FP3>())
        .def("getFields", &pyPSTDGrid::getFields)
        .def("getJ", &pyPSTDGrid::getJ)
        .def("getE", &pyPSTDGrid::getE)
        .def("getB", &pyPSTDGrid::getB)
        .def("setJ", &pyPSTDGrid::setJ)
        .def("setE", &pyPSTDGrid::setE)
        .def("setB", &pyPSTDGrid::setB)
		.def("setJ", &pyPSTDGrid::pySetJ)
		.def("setE", &pyPSTDGrid::pySetE)
		.def("setB", &pyPSTDGrid::pySetB)
        .def("setJ", &pyPSTDGrid::setJxyz)
        .def("setE", &pyPSTDGrid::setExyz)
        .def("setB", &pyPSTDGrid::setBxyz)
		.def("setJ", &pyPSTDGrid::pySetJxyz)
		.def("setE", &pyPSTDGrid::pySetExyz)
		.def("setB", &pyPSTDGrid::pySetBxyz)
		.def("setJt", &pyPSTDGrid::setJxyzt)
		.def("setEt", &pyPSTDGrid::setExyzt)
		.def("setBt", &pyPSTDGrid::setBxyzt)
		.def("analytical", &pyPSTDGrid::setAnalytical)
		.def("setTime", &pyPSTDGrid::setTime)
        ;

    py::class_<PSTD>(object, "PSTD")
        .def(py::init<pyPSTDGrid*>())
        .def("setPML", &PSTD::setPML)
        .def("updateFields", &PSTD::updateFields)
        .def("setTimeStep", &PSTD::setTimeStep)
        ;

    py::class_<pyPSATDGrid>(object, "PSATDGrid")
        .def(py::init<FP3, FP, FP3, FP3>())
        .def("getFields", &pyPSATDGrid::getFields)
        .def("getJ", &pyPSATDGrid::getJ)
        .def("getE", &pyPSATDGrid::getE)
        .def("getB", &pyPSATDGrid::getB)
        .def("setJ", &pyPSATDGrid::setJ)
        .def("setE", &pyPSATDGrid::setE)
        .def("setB", &pyPSATDGrid::setB)
		.def("setJ", &pyPSATDGrid::pySetJ)
		.def("setE", &pyPSATDGrid::pySetE)
		.def("setB", &pyPSATDGrid::pySetB)
        .def("setJ", &pyPSATDGrid::setJxyz)
        .def("setE", &pyPSATDGrid::setExyz)
        .def("setB", &pyPSATDGrid::setBxyz)
		.def("setJ", &pyPSATDGrid::pySetJxyz)
		.def("setE", &pyPSATDGrid::pySetExyz)
		.def("setB", &pyPSATDGrid::pySetBxyz)
		.def("setJt", &pyPSATDGrid::setJxyzt)
		.def("setEt", &pyPSATDGrid::setExyzt)
		.def("setBt", &pyPSATDGrid::setBxyzt)
		.def("analytical", &pyPSATDGrid::setAnalytical)
		.def("setTime", &pyPSATDGrid::setTime)
        ;

    py::class_<PSATD>(object, "PSATD")
        .def(py::init<pyPSATDGrid*>())
        .def("setPML", &PSATD::setPML)
        .def("updateFields", &PSATD::updateFields)
        .def("setTimeStep", &PSATD::setTimeStep)
        ;

    py::class_<PSATDTimeStraggered>(object, "PSATDTimeStraggered")
        .def(py::init<pyPSATDGrid*>())
        .def("setPML", &PSATDTimeStraggered::setPML)
        .def("updateFields", &PSATDTimeStraggered::updateFields)
        .def("setTimeStep", &PSATDTimeStraggered::setTimeStep)
        ;

	py::class_<ScalarQED_AEG_only_electron>(object, "QED")
		.def(py::init<>())
		.def("processParticles", &ScalarQED_AEG_only_electron::processParticles)
		.def("processParticlesNIter", &ScalarQED_AEG_only_electron::processParticlesNIter)
		;
}
