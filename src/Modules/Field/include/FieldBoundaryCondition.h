#pragma once
#include "Enums.h"
#include "Vectors.h"

namespace pfc
{

    template<class TGrid>
    class FieldBoundaryCondition
    {
    public:

        FieldBoundaryCondition(TGrid* grid, CoordinateEnum axis) :
            grid(grid), axis(axis) {
        }

        // polymorfic class
        virtual ~FieldBoundaryCondition() {}
        virtual FieldBoundaryCondition<TGrid>* createInstance(
            TGrid* grid, CoordinateEnum axis) = 0;

        virtual void generateB(FP time) = 0;
        virtual void generateE(FP time) = 0;

        TGrid* grid = nullptr;
        CoordinateEnum axis;
    };
}
