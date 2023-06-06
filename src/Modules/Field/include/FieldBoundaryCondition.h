#pragma once
#include <array>

#include "Grid.h"
#include "FieldSolver.h"

namespace pfc
{

    template<GridTypes gridTypes>
    class FieldBoundaryCondition
    {
    public:

        FieldBoundaryCondition(CoordinateEnum axis, FieldSolver<gridTypes>* fieldSolver) :
            fieldSolver(fieldSolver), axis(axis) {
        }

        FieldBoundaryCondition(const FieldBoundaryCondition<gridTypes>& gen) :
            fieldSolver(gen.fieldSolver), axis(gen.axis) {
        }

        virtual ~FieldBoundaryCondition() {}

        virtual void generateB() = 0;
        virtual void generateE() = 0;

        virtual FieldBoundaryCondition<gridTypes>* createInstance(
            FieldSolver<gridTypes>* fieldSolver = nullptr) = 0;

        FieldSolver<gridTypes>* fieldSolver;
        CoordinateEnum axis;
    };
}
