#pragma once
#include "FieldBoundaryCondition.h"

namespace pfc
{
    template <class TGrid>
    class PeriodicalBoundaryConditionSpectral : public FieldBoundaryCondition<TGrid> {
    public:

        // periodical boundaries for spectral solvers are default because of FFT
        // this class does nothing

        PeriodicalBoundaryConditionSpectral(TGrid* grid, CoordinateEnum axis,
            Int3 leftBorderIndex, Int3 rightBorderIndex) :
            FieldBoundaryCondition<TGrid>(grid, axis, leftBorderIndex, rightBorderIndex) {}

        void generateB(FP time) override {}
        void generateE(FP time) override {}

        FieldBoundaryCondition<TGrid>* createInstance(
            TGrid* grid, CoordinateEnum axis, Int3 leftBorderIndex, Int3 rightBorderIndex) override {
            return new PeriodicalBoundaryConditionSpectral(grid, axis, leftBorderIndex, rightBorderIndex);
        }
    };
}
