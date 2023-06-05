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

        FieldBoundaryCondition(CoordinateEnum axis, FieldSolver<gridTypes>* fieldSolver = 0) :
            fieldSolver(fieldSolver), axis(axis) {
        }

        FieldBoundaryCondition(const FieldBoundaryCondition<gridTypes>& gen,
            FieldSolver<gridTypes>* fieldSolver = 0) :
            fieldSolver(fieldSolver ? fieldSolver: gen.fieldSolver), axis(gen.axis) {
        }

        virtual ~FieldBoundaryCondition() {}

        virtual void generateB() = 0;
        virtual void generateE() = 0;

        virtual FieldBoundaryCondition<gridTypes>* createInstance(FieldSolver<gridTypes>* fieldSolver) = 0;

        FieldSolver<gridTypes>* fieldSolver;
        CoordinateEnum axis;
    };
}
