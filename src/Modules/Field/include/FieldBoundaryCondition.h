#pragma once
#include "Enums.h"
#include "Vectors.h"
#include "Grid.h"

namespace pfc
{

    template<class TGrid>
    class FieldBoundaryCondition
    {
    public:

        FieldBoundaryCondition(TGrid* grid, CoordinateEnum axis,
            Int3 leftBorderIndex, Int3 rightBorderIndex) :
            grid(grid), axis(axis),
            leftBorderIndex(leftBorderIndex), rightBorderIndex(rightBorderIndex) {}

        // polymorfic class
        virtual ~FieldBoundaryCondition() {}
        virtual FieldBoundaryCondition<TGrid>* createInstance(
            TGrid* grid, CoordinateEnum axis, Int3 leftBorderIndex, Int3 rightBorderIndex) = 0;

        virtual void generateB(FP time) = 0;
        virtual void generateE(FP time) = 0;

        TGrid* grid = nullptr;
        CoordinateEnum axis;
        Int3 leftBorderIndex, rightBorderIndex;
    };
}
