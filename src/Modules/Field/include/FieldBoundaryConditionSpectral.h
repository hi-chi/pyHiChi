#pragma once
#include "FieldBoundaryCondition.h"

namespace pfc
{
    template <GridTypes gridTypes>
    class SpectralPeriodicalBoundaryCondition : public FieldBoundaryCondition<gridTypes> {
    public:

        // periodical boundaries for spectral solvers are default because of FFT
        // this class does nothing

        SpectralPeriodicalBoundaryCondition(SpectralFieldSolver<gridTypes>* fieldSolver = 0,
            bool isXAxisEnabled = true, bool isYAxisEnabled = true, bool isZAxisEnabled = true) :
            FieldBoundaryCondition<gridTypes>(fieldSolver, isXAxisEnabled, isYAxisEnabled, isZAxisEnabled) {
        }

        SpectralPeriodicalBoundaryCondition(const SpectralPeriodicalBoundaryCondition& gen,
            SpectralFieldSolver<gridTypes>* fieldSolver = 0) :
            FieldBoundaryCondition<gridTypes>(gen, fieldSolver) {
        }

        void generateB() override {};
        void generateE() override {};

        FieldBoundaryCondition<gridTypes>* createInstance(FieldSolver<gridTypes>* fieldSolver) override {
            return new SpectralPeriodicalBoundaryCondition(*this, static_cast<SpectralFieldSolver<gridTypes>*>(fieldSolver));
        }
    };

    typedef SpectralPeriodicalBoundaryCondition<GridTypes::PSTDGridType> PeriodicalBoundaryConditionPstd;
    typedef SpectralPeriodicalBoundaryCondition<GridTypes::PSATDGridType> PeriodicalBoundaryConditionPsatd;
    typedef SpectralPeriodicalBoundaryCondition<GridTypes::PSATDTimeStraggeredGridType> PeriodicalBoundaryConditionPsatdTimeStraggered;
}