#pragma once
#include "FieldBoundaryCondition.h"

namespace pfc
{
    template <class TGrid>
    class PeriodicalBoundaryConditionSpectral : public FieldBoundaryCondition<TGrid> {
    public:

        // periodical boundaries for spectral solvers are default because of FFT
        // this class does nothing

        PeriodicalBoundaryConditionSpectral(TGrid* grid, CoordinateEnum axis) :
            FieldBoundaryCondition<TGrid>(grid, axis) {
        }

        void generateB(FP time) override {}
        void generateE(FP time) override {}

        FieldBoundaryCondition<TGrid>* createInstance(
            TGrid* grid, CoordinateEnum axis) override {
            return new PeriodicalBoundaryConditionSpectral(grid, axis);
        }
    };
}
