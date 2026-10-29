#pragma once

#include "pyField.h"
#include "CurrentDeposition.h"
#include "CurrentBoundaries.h"
#include "ParticleGenerator.h"
#include "ParticleBoundaries.h"

namespace pfc
{
    template<class TFieldSolver>
    class pyModuleInterface {
    public:
        pyModuleInterface(pyField<TFieldSolver>* _field): field(_field) {}

    protected:
        pyField<TFieldSolver>* field;
    };


    template<class TFieldSolver>
    class pyParticleGenerator: public ParticleGenerator,
                               public pyModuleInterface<TFieldSolver> {
    public:
        typedef FP WeightType;
        typedef ParticleTypes TypeIndexType;

        pyParticleGenerator(pyField<TFieldSolver>* _field):
            ParticleGenerator(), pyModuleInterface<TFieldSolver>(_field) {}

        template <class TParticleArray>
        void operator()(TParticleArray* particleArray,
            int64_t particleDensity,
            int64_t initialTemperature,
            WeightType weight = 1.0,
            TypeIndexType typeIndex = ParticleTypes::Electron)
        {
            ParticleGenerator::operator()(
                particleArray, this->field->getGrid(),
                (FP(*)(FP, FP, FP))particleDensity,
                (FP(*)(FP, FP, FP))initialTemperature,
                [](FP init_mx, FP init_my, FP init_mz) -> 
                    FP3 {return FP3(0.0, 0.0, 0.0); },
                weight, typeIndex);
        }
    };

    template<class TGrid, class TFieldSolver>
    class pyCurrentDepositionCIC: public CurrentDepositionCIC<TGrid>,
                                  public pyModuleInterface<TFieldSolver> {
    public:
        pyCurrentDepositionCIC(pyField<TFieldSolver>* _field, double _dt):
            CurrentDepositionCIC<TGrid>(_dt), pyModuleInterface<TFieldSolver>(_field) {}

        template<class TParticleArray>
        void operator()(TParticleArray* particleArray) {
            CurrentDepositionCIC<TGrid>::operator()(this->field->getGrid(), particleArray);
        }

        FP3 getJ(const Int3& idx) {
            return this->field->getJ(idx);
        }
    };

    template<class TFieldSolver>
    class pyPeriodicalParticleBoundaryConditions : public PeriodicalParticleBoundaryConditions,
                                                   public pyModuleInterface<TFieldSolver> {
    public:
        pyPeriodicalParticleBoundaryConditions(pyField<TFieldSolver>* _field):
            PeriodicalParticleBoundaryConditions(), pyModuleInterface<TFieldSolver>(_field) {}

        template<class TParticleArray>
        void update(TParticleArray* particleArray) {
            PeriodicalParticleBoundaryConditions::updateParticlePosition(
                this->field->getGrid(), particleArray);
        }
    };

    template<class TFieldSolver>
    class pyInterpolation : public pyModuleInterface<TFieldSolver> {
    public:
        pyInterpolation(pyField<TFieldSolver>* _field) : pyModuleInterface<TFieldSolver>(_field) {}

        FP getExCIC(const FP3& coords) {
            return this->field->getGrid()->getExCIC(coords);
        }
        FP getEyCIC(const FP3& coords) {
            return this->field->getGrid()->getEyCIC(coords);
        }
        FP getEzCIC(const FP3& coords) {
            return this->field->getGrid()->getEzCIC(coords);
        }
        FP getBxCIC(const FP3& coords) {
            return this->field->getGrid()->getBxCIC(coords);
        }
        FP getByCIC(const FP3& coords) {
            return this->field->getGrid()->getByCIC(coords);
        }
        FP getBzCIC(const FP3& coords) {
            return this->field->getGrid()->getBzCIC(coords);
        }
        FP getJxCIC(const FP3& coords) {
            return this->field->getGrid()->getJxCIC(coords);
        }
        FP getJyCIC(const FP3& coords) {
            return this->field->getGrid()->getJyCIC(coords);
        }
        FP getJzCIC(const FP3& coords) {
            return this->field->getGrid()->getJzCIC(coords);
        }
    };


    template<class TFieldSolver>
    class pyPeriodicalCurrentBC: public PeriodicalCurrentBoundaryConditionFdtd,
                                 public pyModuleInterface<TFieldSolver> {
    public:
        pyPeriodicalCurrentBC(pyField<TFieldSolver>* _field) : 
            PeriodicalCurrentBoundaryConditionFdtd(_field->getGrid()),
            pyModuleInterface<TFieldSolver>(_field) {}
        void update() {
            this->field->getFieldSolver()->updateDomainBorders();
            PeriodicalCurrentBoundaryConditionFdtd::updateCurrentBoundaries();
        }
    };
}
