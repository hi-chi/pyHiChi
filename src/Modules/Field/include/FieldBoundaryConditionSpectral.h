#pragma once
#include "FieldBoundaryCondition.h"

namespace pfc
{
    template <class TGrid>
    class PeriodicalBoundaryConditionSpectral : public FieldBoundaryCondition<TGrid> {
    public:

        PeriodicalBoundaryConditionSpectral(TGrid* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) :
            FieldBoundaryCondition<TGrid>(grid, leftBorderIndex, rightBorderIndex, axis)
        {}

        // constructor for loading
        PeriodicalBoundaryConditionSpectral(TGrid* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex) :
            FieldBoundaryCondition<TGrid>(grid, leftBorderIndex, rightBorderIndex)
        {}

        // periodical boundaries for spectral solvers are default because of FFT
        // this class does nothing
        void generateB(FP time) override {}
        void generateE(FP time) override {}

        FieldBoundaryCondition<TGrid>* createInstance(
            TGrid* grid, Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) override {
            return new PeriodicalBoundaryConditionSpectral(grid, leftBorderIndex, rightBorderIndex, axis);
        }
    };
}
