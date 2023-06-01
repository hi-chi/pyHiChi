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

        FieldBoundaryCondition(FieldSolver<gridTypes>* fieldSolver = 0,
            bool isXAxisEnabled = true, bool isYAxisEnabled = true, bool isZAxisEnabled = true) :
            fieldSolver(fieldSolver),
            enabledAxis({ isXAxisEnabled, isYAxisEnabled, isZAxisEnabled }) {
        }

        // copy constructor, other fieldSolver is possible
        FieldBoundaryCondition(const FieldBoundaryCondition<gridTypes>& gen,
            FieldSolver<gridTypes>* fieldSolver = 0) :
            fieldSolver(gen.fieldSolver), enabledAxis(gen.enabledAxis) {
        }

        virtual ~FieldBoundaryCondition() {}

        virtual void generateB() = 0;
        virtual void generateE() = 0;

        virtual FieldBoundaryCondition<gridTypes>* createInstance(FieldSolver<gridTypes>* fieldSolver) = 0;

        FieldSolver<gridTypes>* fieldSolver;
        std::array<bool, 3> enabledAxis;
    };
}
