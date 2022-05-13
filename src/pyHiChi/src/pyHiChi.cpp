// Python-Interface.cpp : Defines exported functions for the dll-file.
//

#include <algorithm>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <pybind11/operators.h>

#include "pyField.h"
#include "pyFieldMacroses.h"

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


namespace py = pybind11;
using namespace pfc;


std::vector<ParticleType> ParticleInfo::typesVector = { {constants::electronMass, constants::electronCharge},//electron
                                    {constants::electronMass, -constants::electronCharge},//positron
                                    {constants::protonMass, -constants::electronCharge},//proton
                                    {constants::electronMass, 0.0 } };//photon
const ParticleType* ParticleInfo::types = &ParticleInfo::typesVector[0];
short ParticleInfo::numTypes = sizeParticleTypes;

template <class QED, class Field, class Grid>
void processParticles(QED* self, Ensemble3d* particles,
    Field* field, FP timeStep, FP startTime, int N)
{
    for (int i = 0; i < N; i++)
    {
        field->setTime(startTime + i * timeStep);
        self->processParticles(particles,
            field->getField()->getGrid(), timeStep);
    }
}


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

    py::enum_<CoordinateEnum>(object, "Axis")
        .value("X", CoordinateEnum::x)
        .value("Y", CoordinateEnum::y)
        .value("Z", CoordinateEnum::z)
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

    py::class_<ValueField>(object, "FieldValue")
        .def(py::init<FP3, FP3>(), py::arg("E"), py::arg("B"))
        .def(py::init<FP, FP, FP, FP, FP, FP>(),
            py::arg("Ex"), py::arg("Ey"), py::arg("Ez"),
            py::arg("Bx"), py::arg("By"), py::arg("Bz"))
        .def("get_E", &ValueField::getE)
        .def("set_E", &ValueField::setE, py::arg("E"))
        .def("get_B", &ValueField::getB)
        .def("set_B", &ValueField::setB, py::arg("B"))
        .def_readwrite("E", &ValueField::E)
        .def_readwrite("B", &ValueField::B)
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

    // -------------------------- QED ---------------------------

    py::class_<ScalarQED_AEG_only_electron_Yee>(object, "QED_Yee")
        .def(py::init<>())
        .def("process_particles", &ScalarQED_AEG_only_electron_Yee::processParticles)
        .def("process_particles", &processParticles<ScalarQED_AEG_only_electron_Yee,
            pyYeeField, YeeGrid>)
        ;

    py::class_<ScalarQED_AEG_only_electron_PSTD>(object, "QED_PSTD")
        .def(py::init<>())
        .def("process_particles", &ScalarQED_AEG_only_electron_PSTD::processParticles)
        .def("process_particles", &processParticles<ScalarQED_AEG_only_electron_PSTD,
            pyPSTDField, PSTDGrid>)
        ;

    py::class_<ScalarQED_AEG_only_electron_PSATD>(object, "QED_PSATD")
        .def(py::init<>())
        .def("process_particles", &ScalarQED_AEG_only_electron_PSATD::processParticles)
        .def("process_particles", &processParticles<ScalarQED_AEG_only_electron_PSATD,
            pyPSATDField, PSATDGrid>)
        ;

    py::class_<ScalarQED_AEG_only_electron_Analytical>(object, "QED_Analytical")
        .def(py::init<>())
        .def("process_particles", &ScalarQED_AEG_only_electron_Analytical::processParticles)
        .def("process_particles", &processParticles<ScalarQED_AEG_only_electron_Analytical,
            pyAnalyticalField, AnalyticalField>)
        ;

    // ------------------- thinnings -------------------
   
    py::enum_<Thinning<ParticleArray3d>::Features>(object, "Conserve")
        .value("momentum", Thinning<ParticleArray3d>::Features::Momentum)
        .value("position", Thinning<ParticleArray3d>::Features::Position)
        .value("energy", Thinning<ParticleArray3d>::Features::Energy)
        .value("dis_momentum", Thinning<ParticleArray3d>::Features::Dispersion_Momentum)
        .value("dis_position", Thinning<ParticleArray3d>::Features::Dispersion_Position)
        .value("dis_energy", Thinning<ParticleArray3d>::Features::Dispersion_Energy)
        .export_values()
        ; 

    py::class_<Thinning<ParticleArray3d>>(object, "Thinout")
        .def(py::init<>())
        .def("simple_thinning", &Thinning<ParticleArray3d>::simple)
        .def("leveling_thinning", &Thinning<ParticleArray3d>::leveling)
        .def("number_conservative_thinning", &Thinning<ParticleArray3d>::numberConservative)
        .def("energy_conservative_thinning", &Thinning<ParticleArray3d>::energyConservative)
        .def("conservative_thinning", &Thinning<ParticleArray3d>::thinningConservative)
        .def("k_means_mergining", &Merging<ParticleArray3d>::merge_with_kmeans)
        ; 

    // ------------------- mappings -------------------

    py::class_<Mapping, std::shared_ptr<Mapping>> pyMapping(object, "Mapping");

    py::class_<IdentityMapping, std::shared_ptr<IdentityMapping>>(object, "IdentityMapping", pyMapping)
        .def(py::init<const FP3&, const FP3&>(), py::arg("a"), py::arg("b"))
        .def("get_direct_coords", &IdentityMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &IdentityMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<PeriodicalMapping, std::shared_ptr<PeriodicalMapping>>(object, "PeriodicalMapping", pyMapping)
        .def(py::init<CoordinateEnum, FP, FP>(), py::arg("axis"), py::arg("c_min"), py::arg("c_max"))
        .def("get_direct_coords", &PeriodicalMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &PeriodicalMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;
    
    py::class_<RotationMapping, std::shared_ptr<RotationMapping>>(object, "RotationMapping", pyMapping)
        .def(py::init<CoordinateEnum, FP>(), py::arg("axis"), py::arg("angle"))
        .def(py::init<CoordinateEnum, FP, CoordinateEnum, FP>(), py::arg("axis"), py::arg("angle"),
            py::arg("prop_dir"), py::arg("pol_angle"))
        .def("get_direct_coords", &RotationMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &RotationMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<ScaleMapping, std::shared_ptr<ScaleMapping>>(object, "ScaleMapping", pyMapping)
        .def(py::init<CoordinateEnum, FP>(), py::arg("axis"), py::arg("scale"))
        .def("get_direct_coords", &ScaleMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &ScaleMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<ShiftMapping, std::shared_ptr<ShiftMapping>>(object, "ShiftMapping", pyMapping)
        .def(py::init<FP3>(), py::arg("shift"))
        .def("get_direct_coords", &ShiftMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &ShiftMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        ;

    py::class_<TightFocusingMapping, std::shared_ptr<TightFocusingMapping>>(object, "TightFocusingMapping", pyMapping)
        .def(py::init<FP, FP, FP>(), py::arg("R0"), py::arg("L"), py::arg("D"))
        .def(py::init<FP, FP, FP, CoordinateEnum>(), py::arg("R0"), py::arg("L"),
            py::arg("D"), py::arg("axis"))
        .def("get_direct_coords", &TightFocusingMapping::getDirectCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_inverse_coords", &TightFocusingMapping::getInverseCoords, py::arg("coords"),
            py::arg("time") = 0.0, py::arg("status") = 0)
        .def("get_min_coord", &TightFocusingMapping::getMinCoord)
        .def("get_max_coord", &TightFocusingMapping::getMaxCoord)
        .def("if_perform_inverse_mapping", &TightFocusingMapping::setIfCut)
        ;

    // ------------------- py scalar field -------------------

    py::class_<pyScalarField, std::shared_ptr<pyScalarField>>(
        object, "ScalarField", py::buffer_protocol())
        .def_buffer([](pyScalarField& sf) -> py::buffer_info {  // cast to np.array
                return py::buffer_info(
                    sf.getData(),                               // Pointer to buffer
                    sizeof(FP),                                 // Size of one scalar
                    py::format_descriptor<FP>::format(),        // Python struct-style format descriptor
                    3,                                          // Number of dimensions
                    { sf.getSize().x, sf.getSize().y, sf.getSize().z },  // Buffer dimensions
                    {                                                    // Strides (in bytes) for each index
                        sizeof(FP) * sf.getSize().y * sf.getSize().z,    
                        sizeof(FP) * sf.getSize().z,
                        sizeof(FP)
                    }
                );
            })
        .def("get_size", &pyScalarField::getSize)
        .def("get", static_cast<FP(pyScalarField::*)(int, int, int) const>(&pyScalarField::get))
        ;

    // ------------------- py fields -------------------

    // abstract class
    py::class_<pyFieldBase, std::shared_ptr<pyFieldBase>> pyClassFieldBase(object, "FieldBase");
    pyClassFieldBase.def("get_fields", &pyFieldBase::getFields)
        .def("get_E", &pyFieldBase::getE, py::arg("coords"))
        .def("get_B", &pyFieldBase::getB, py::arg("coords"))
        .def("get_J", &pyFieldBase::getJ, py::arg("coords"))
        .def("get_Ex", &pyFieldBase::getEx, py::arg("coords"))
        .def("get_Ey", &pyFieldBase::getEy, py::arg("coords"))
        .def("get_Ez", &pyFieldBase::getEz, py::arg("coords"))
        .def("get_Bx", &pyFieldBase::getBx, py::arg("coords"))
        .def("get_By", &pyFieldBase::getBy, py::arg("coords"))
        .def("get_Bz", &pyFieldBase::getBz, py::arg("coords"))
        .def("get_Jx", &pyFieldBase::getJx, py::arg("coords"))
        .def("get_Jy", &pyFieldBase::getJy, py::arg("coords"))
        .def("get_Jz", &pyFieldBase::getJz, py::arg("coords"))
        .def("get_E", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getE(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_B", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getB(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_J", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getJ(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_Ex", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getEx(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_Ey", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getEy(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_Ez", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getEz(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_Bx", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getBx(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_By", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getBy(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_Bz", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getBz(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_Jx", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getJx(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_Jy", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getJy(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_Jz", [](std::shared_ptr<pyFieldBase> self, FP x, FP y, FP z) {
                return self->getJz(FP3(x, y, z));
            }, py::arg("x"), py::arg("y"), py::arg("z"))
        
        .def("update_fields", &pyFieldBase::updateFields)
        .def("advance", &pyFieldBase::advance, py::arg("time_step"))

        SET_ALL_PYFIELDBASE_SECTION_METHODS()
        ;

    // sum and mul fields

    py::class_<pySumField, std::shared_ptr<pySumField>>(
        object, "SumField", pyClassFieldBase)
        SET_SUM_AND_MAP_FIELD_METHODS(pySumField)
        ;

    py::class_<pyMulField, std::shared_ptr<pyMulField>>(
        object, "MulField", pyClassFieldBase)
        SET_SUM_AND_MAP_FIELD_METHODS(pyMulField)
        ;

    // simple fields

    py::class_<pyAnalyticalField, std::shared_ptr<pyAnalyticalField>>(
        object, "AnalyticalField", pyClassFieldBase)
        SET_SUM_AND_MAP_FIELD_METHODS(pyAnalyticalField)
        SET_COMMON_FIELD_METHODS(pyAnalyticalField)
        .def(py::init<FP>(), py::arg("time_step"))
        .def("set_E", &pyAnalyticalField::setExyz,
            py::arg("Ex"), py::arg("Ey"), py::arg("Ez"))
        .def("set_B", &pyAnalyticalField::setBxyz,
            py::arg("Bx"), py::arg("By"), py::arg("Bz"))
        .def("set_J", &pyAnalyticalField::setJxyz,
            py::arg("Jx"), py::arg("Jy"), py::arg("Jz"))
        .def("get_E", &pyAnalyticalField::getEt,
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("t"))
        .def("get_B", &pyAnalyticalField::getBt,
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("t"))
        .def("get_J", &pyAnalyticalField::getJt,
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("t"))
        ;

    py::class_<pyYeeField, std::shared_ptr<pyYeeField>>(
        object, "YeeField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>(), py::arg("grid_size"), py::arg("min_coords"),
            py::arg("spatial_steps"), py::arg("time_step"))
        SET_COMPUTATIONAL_GRID_METHODS(pyYeeField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyYeeField)
        SET_SCALAR_FIELD_METHODS(pyYeeField)
        SET_COMMON_FIELD_METHODS(pyYeeField)
        .def("set_PML", &pyYeeField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("set_periodical_BC", &pyYeeField::setPeriodicalFieldGenerator)
        .def("zoom", &pyYeeField::zoom, py::arg("min_coord"), \
            py::arg("zoomed_grid_size"), py::arg("zoomed_grid_step"))
        ;

    py::class_<pyPSTDField, std::shared_ptr<pyPSTDField>>(
        object, "PSTDField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>(), py::arg("grid_size"), py::arg("min_coords"),
            py::arg("spatial_steps"), py::arg("time_step"))
        SET_COMPUTATIONAL_GRID_METHODS(pyPSTDField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyPSTDField)
        SET_SCALAR_FIELD_METHODS(pyPSTDField)
        SET_COMMON_FIELD_METHODS(pyPSTDField)
        .def("set_PML", &pyPSTDField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("set", &pyPSTDField::setEMField, py::arg("func"))
        .def("set", &pyPSTDField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyPSTDField::applyFunction, py::arg("func"))
        .def("apply_function", &pyPSTDField::pyApplyFunction, py::arg("func"))
        .def("zoom", &pyPSTDField::zoom, py::arg("min_coord"), \
            py::arg("zoomed_grid_size"), py::arg("zoomed_grid_step"))
        ;

    py::class_<pyPSATDField, std::shared_ptr<pyPSATDField>>(
        object, "PSATDField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>(), py::arg("grid_size"), py::arg("min_coords"),
            py::arg("spatial_steps"), py::arg("time_step"))
        SET_COMPUTATIONAL_GRID_METHODS(pyPSATDField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyPSATDField)
        SET_SCALAR_FIELD_METHODS(pyPSATDField)
        SET_COMMON_FIELD_METHODS(pyPSATDField)
        .def("set_PML", &pyPSATDField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("convert_fields_poisson_equation", &pyPSATDField::convertFieldsPoissonEquation)
        .def("set", &pyPSATDField::setEMField, py::arg("func"))
        .def("set", &pyPSATDField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyPSATDField::applyFunction, py::arg("func"))
        .def("apply_function", &pyPSATDField::pyApplyFunction, py::arg("func"))
        .def("zoom", &pyPSATDField::zoom, py::arg("min_coord"), \
            py::arg("zoomed_grid_size"), py::arg("zoomed_grid_step"))
        ;

    py::class_<pyPSATDPoissonField, std::shared_ptr<pyPSATDPoissonField>>(
        object, "PSATDPoissonField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>(), py::arg("grid_size"), py::arg("min_coords"),
            py::arg("spatial_steps"), py::arg("time_step"))
        SET_COMPUTATIONAL_GRID_METHODS(pyPSATDPoissonField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyPSATDPoissonField)
        SET_SCALAR_FIELD_METHODS(pyPSATDPoissonField)
        SET_COMMON_FIELD_METHODS(pyPSATDPoissonField)
        .def("set_PML", &pyPSATDPoissonField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("convert_fields_poisson_equation", &pyPSATDPoissonField::convertFieldsPoissonEquation)
        .def("set", &pyPSATDPoissonField::setEMField, py::arg("func"))
        .def("set", &pyPSATDPoissonField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyPSATDPoissonField::applyFunction, py::arg("func"))
        .def("apply_function", &pyPSATDPoissonField::pyApplyFunction, py::arg("func"))
        .def("zoom", &pyPSATDPoissonField::zoom, py::arg("min_coord"), \
            py::arg("zoomed_grid_size"), py::arg("zoomed_grid_step"))
        ;

    py::class_<pyPSATDTimeStraggeredField, std::shared_ptr<pyPSATDTimeStraggeredField>>(
        object, "PSATDSField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>(), py::arg("grid_size"), py::arg("min_coords"),
            py::arg("spatial_steps"), py::arg("time_step"))
        SET_COMPUTATIONAL_GRID_METHODS(pyPSATDTimeStraggeredField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyPSATDTimeStraggeredField)
        SET_SCALAR_FIELD_METHODS(pyPSATDTimeStraggeredField)
        SET_COMMON_FIELD_METHODS(pyPSATDTimeStraggeredField)
        .def("set_PML", &pyPSATDTimeStraggeredField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("convert_fields_poisson_equation", &pyPSATDTimeStraggeredField::convertFieldsPoissonEquation)
        .def("set", &pyPSATDTimeStraggeredField::setEMField, py::arg("func"))
        .def("set", &pyPSATDTimeStraggeredField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyPSATDTimeStraggeredField::applyFunction, py::arg("func"))
        .def("apply_function", &pyPSATDTimeStraggeredField::pyApplyFunction, py::arg("func"))
        .def("zoom", &pyPSATDTimeStraggeredField::zoom, py::arg("min_coord"), \
            py::arg("zoomed_grid_size"), py::arg("zoomed_grid_step"))
        ;

    py::class_<pyPSATDTimeStraggeredPoissonField, std::shared_ptr<pyPSATDTimeStraggeredPoissonField>>(
        object, "PSATDSPoissonField", pyClassFieldBase)
        .def(py::init<FP3, FP3, FP3, FP>(), py::arg("grid_size"), py::arg("min_coords"),
            py::arg("spatial_steps"), py::arg("time_step"))
        SET_COMPUTATIONAL_GRID_METHODS(pyPSATDTimeStraggeredPoissonField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyPSATDTimeStraggeredPoissonField)
        SET_SCALAR_FIELD_METHODS(pyPSATDTimeStraggeredPoissonField)
        SET_COMMON_FIELD_METHODS(pyPSATDTimeStraggeredPoissonField)
        .def("set_PML", &pyPSATDTimeStraggeredPoissonField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("convert_fields_poisson_equation", &pyPSATDTimeStraggeredPoissonField::convertFieldsPoissonEquation)
        .def("set", &pyPSATDTimeStraggeredPoissonField::setEMField, py::arg("func"))
        .def("set", &pyPSATDTimeStraggeredPoissonField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyPSATDTimeStraggeredPoissonField::applyFunction, py::arg("func"))
        .def("apply_function", &pyPSATDTimeStraggeredPoissonField::pyApplyFunction, py::arg("func"))
        .def("zoom", &pyPSATDTimeStraggeredPoissonField::zoom, py::arg("min_coord"), \
            py::arg("zoomed_grid_size"), py::arg("zoomed_grid_step"))
        ;

    // mapped fields

    py::class_<pyMappedAnalyticalField, std::shared_ptr<pyMappedAnalyticalField>>(
        object, "MappedAnalyticalField", pyClassFieldBase)
        SET_SUM_AND_MAP_FIELD_METHODS(pyMappedAnalyticalField)
        SET_COMMON_FIELD_METHODS(pyMappedAnalyticalField)
        .def("set_E", &pyMappedAnalyticalField::setExyz,
            py::arg("Ex"), py::arg("Ey"), py::arg("Ez"))
        .def("set_B", &pyMappedAnalyticalField::setBxyz,
            py::arg("Bx"), py::arg("By"), py::arg("Bz"))
        .def("set_J", &pyMappedAnalyticalField::setJxyz,
            py::arg("Jx"), py::arg("Jy"), py::arg("Jz"))
        .def("get_E", &pyMappedAnalyticalField::getEt,
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("t"))
        .def("get_B", &pyMappedAnalyticalField::getBt,
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("t"))
        .def("get_J", &pyMappedAnalyticalField::getJt,
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("t"))
        ;

    py::class_<pyMappedYeeField, std::shared_ptr<pyMappedYeeField>>(
        object, "MappedYeeField", pyClassFieldBase)
        SET_COMPUTATIONAL_GRID_METHODS(pyMappedYeeField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyMappedYeeField)
        SET_COMMON_FIELD_METHODS(pyMappedYeeField)
        .def("set_PML", &pyMappedYeeField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("set_periodical_BC", &pyMappedYeeField::setPeriodicalFieldGenerator)
        ;

    py::class_<pyMappedPSTDField, std::shared_ptr<pyMappedPSTDField>>(
        object, "MappedPSTDField", pyClassFieldBase)
        SET_COMPUTATIONAL_GRID_METHODS(pyMappedPSTDField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyMappedPSTDField)
        SET_COMMON_FIELD_METHODS(pyMappedPSTDField)
        .def("set_PML", &pyMappedPSTDField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("set", &pyMappedPSTDField::setEMField, py::arg("func"))
        .def("set", &pyMappedPSTDField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyMappedPSTDField::applyFunction, py::arg("func"))
        .def("apply_function", &pyMappedPSTDField::pyApplyFunction, py::arg("func"))
        ;

    py::class_<pyMappedPSATDField, std::shared_ptr<pyMappedPSATDField>>(
        object, "MappedPSATDField", pyClassFieldBase)
        SET_COMPUTATIONAL_GRID_METHODS(pyMappedPSATDField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyMappedPSATDField)
        SET_COMMON_FIELD_METHODS(pyMappedPSATDField)
        .def("set_PML", &pyMappedPSATDField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("convert_fields_poisson_equation", &pyMappedPSATDField::convertFieldsPoissonEquation)
        .def("set", &pyMappedPSATDField::setEMField, py::arg("func"))
        .def("set", &pyMappedPSATDField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyMappedPSATDField::applyFunction, py::arg("func"))
        .def("apply_function", &pyMappedPSATDField::pyApplyFunction, py::arg("func"))
        ;

    py::class_<pyMappedPSATDPoissonField, std::shared_ptr<pyMappedPSATDPoissonField>>(
        object, "MappedPSATDPoissonField", pyClassFieldBase)
        SET_COMPUTATIONAL_GRID_METHODS(pyMappedPSATDPoissonField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyMappedPSATDPoissonField)
        SET_COMMON_FIELD_METHODS(pyMappedPSATDPoissonField)
        .def("set_PML", &pyMappedPSATDPoissonField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("convert_fields_poisson_equation", &pyMappedPSATDPoissonField::convertFieldsPoissonEquation)
        .def("set", &pyMappedPSATDPoissonField::setEMField, py::arg("func"))
        .def("set", &pyMappedPSATDPoissonField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyMappedPSATDPoissonField::applyFunction, py::arg("func"))
        .def("apply_function", &pyMappedPSATDPoissonField::pyApplyFunction, py::arg("func"))
        ;

    py::class_<pyMappedPSATDTimeStraggeredField, std::shared_ptr<pyMappedPSATDTimeStraggeredField>>(
        object, "MappedPSATDSField", pyClassFieldBase)
        SET_COMPUTATIONAL_GRID_METHODS(pyMappedPSATDTimeStraggeredField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyMappedPSATDTimeStraggeredField)
        SET_COMMON_FIELD_METHODS(pyMappedPSATDTimeStraggeredField)
        .def("set_PML", &pyMappedPSATDTimeStraggeredField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("convert_fields_poisson_equation", &pyMappedPSATDTimeStraggeredField::convertFieldsPoissonEquation)
        .def("set", &pyMappedPSATDTimeStraggeredField::setEMField, py::arg("func"))
        .def("set", &pyMappedPSATDTimeStraggeredField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyMappedPSATDTimeStraggeredField::applyFunction, py::arg("func"))
        .def("apply_function", &pyMappedPSATDTimeStraggeredField::pyApplyFunction, py::arg("func"))
        ;

    py::class_<pyMappedPSATDTimeStraggeredPoissonField, std::shared_ptr<pyMappedPSATDTimeStraggeredPoissonField>>(
        object, "MappedPSATDSPoissonField", pyClassFieldBase)
        SET_COMPUTATIONAL_GRID_METHODS(pyMappedPSATDTimeStraggeredPoissonField)
        SET_SUM_AND_MAP_FIELD_METHODS(pyMappedPSATDTimeStraggeredPoissonField)
        SET_COMMON_FIELD_METHODS(pyMappedPSATDTimeStraggeredPoissonField)
        .def("set_PML", &pyMappedPSATDTimeStraggeredPoissonField::setPML,
            py::arg("pml_size_x"), py::arg("pml_size_y"), py::arg("pml_size_z"))
        .def("convert_fields_poisson_equation", &pyMappedPSATDTimeStraggeredPoissonField::convertFieldsPoissonEquation)
        .def("set", &pyMappedPSATDTimeStraggeredPoissonField::setEMField, py::arg("func"))
        .def("set", &pyMappedPSATDTimeStraggeredPoissonField::pySetEMField, py::arg("func"))
        .def("apply_function", &pyMappedPSATDTimeStraggeredPoissonField::applyFunction, py::arg("func"))
        .def("apply_function", &pyMappedPSATDTimeStraggeredPoissonField::pyApplyFunction, py::arg("func"))
        ;

    // ------------------- field configurations -------------------

    py::class_<NullField>(object, "NullField")
        .def(py::init<>())
        .def("get_E", &NullField::getE)
        .def("get_B", &NullField::getB)
        ;

    py::class_<TightFocusingField>(object, "TightFocusingField")
        .def(py::init<FP, FP, FP, FP, FP, FP>(), 
            py::arg("f_number"), py::arg("R0"), py::arg("wavelength"), py::arg("pulselength"),
            py::arg("total_power"), py::arg("edge_smoothing_angle"))
        .def(py::init<FP, FP, FP, FP, FP, FP, FP3>(),
            py::arg("f_number"), py::arg("R0"), py::arg("wavelength"), py::arg("pulselength"),
            py::arg("total_power"), py::arg("edge_smoothing_angle"), py::arg("polarisation"))
        .def("get_E", &TightFocusingField::getE)
        .def("get_B", &TightFocusingField::getB)
        ;

}
