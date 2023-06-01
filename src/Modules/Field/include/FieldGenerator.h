#pragma once
#include "Grid.h"
#include "Vectors.h"

#include <functional>

namespace pfc
{
    template<GridTypes gridTypes>
    class FieldSolver;

    template<GridTypes gridTypes>
    class FieldGenerator
    {
    public:

        using FunctionType = std::function<FP(FP, FP, FP, FP)>;

        FieldGenerator(FieldSolver<gridTypes>* fieldSolver = 0);

        virtual ~FieldGenerator() {}

        // copy constructor, other fieldSolver is possible
        FieldGenerator(const FieldGenerator& gen, FieldSolver<gridTypes>* fieldSolver = 0);

        virtual void generateB() = 0;
        virtual void generateE() = 0;

        FieldSolver<gridTypes>* fieldSolver;

    protected:

        // major index is index of edge, minor index is index of component
        FunctionType eLeft[3][3];
        FunctionType eRight[3][3];
        FunctionType bLeft[3][3];
        FunctionType bRight[3][3];
        FP3 leftCoeff;
        FP3 rightCoeff;
    };

    template<GridTypes gridTypes>
    FieldGenerator<gridTypes>::FieldGenerator(FieldSolver<gridTypes>* _fieldSolver)
    {
        fieldSolver = _fieldSolver;
        leftCoeff = FP3(0, 0, 0);
        rightCoeff = FP3(0, 0, 0);
        for (int d = 0; d < 3; ++d)
        {
            //in seq if fieldGeneration in area
            leftCoeff[d] = 1;
            rightCoeff[d] = 1;
        }
    }

    template<GridTypes gridTypes>
    inline FieldGenerator<gridTypes>::FieldGenerator(const FieldGenerator & gen,
        FieldSolver<gridTypes>* fieldSolver)
    {
        if (fieldSolver)
            this->fieldSolver = fieldSolver;
        else this->fieldSolver = gen.fieldSolver;

        leftCoeff = gen.leftCoeff;
        rightCoeff = gen.rightCoeff;

        for (int f = 0; f < 3; ++f)
            for (int d = 0; d < 3; ++d) {
                eLeft[f][d] = gen.eLeft[f][d];
                eRight[f][d] = gen.eRight[f][d];
                bLeft[f][d] = gen.bLeft[f][d];
                bRight[f][d] = gen.bRight[f][d];
            }
    }
}
