#pragma once
#include "FieldGenerator.h"
#include "FieldSolver.h"

namespace pfc
{
    template<GridTypes gridTypes>
    class FieldGeneratorSpectral : public FieldGenerator<gridTypes>
    {
    public:

        // empty class

        FieldGeneratorSpectral(FieldSolver<gridTypes>* fieldSolver = 0) :
            FieldGenerator<gridTypes>(fieldSolver) {}

        // copy constructor, other fieldSolver is possible
        FieldGeneratorSpectral(const FieldGeneratorSpectral& gen,
            FieldSolver<gridTypes>* fieldSolver = 0) :
            FieldGenerator<gridTypes>(gen, fieldSolver) {}

        void generateB() override {}
        void generateE() override {}
    };
}
