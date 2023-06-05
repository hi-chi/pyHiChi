#pragma once
#include "FieldGenerator.h"
#include "FieldSolver.h"

namespace pfc
{
    // temporary empty implementation of field generator for spectral field solvers
    template<GridTypes gridTypes>
    class FieldGeneratorSpectral : public FieldGenerator<gridTypes>
    {
    public:

        FieldGeneratorSpectral(FieldSolver<gridTypes>* fieldSolver,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
            FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1)) :
            FieldGenerator(fieldSolver, leftGenIndex, rightGenIndex,
                bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
                isLeftBorderEnabled, isRightBorderEnabled) {}

        FieldGeneratorSpectral(FieldSolver<gridTypes>* fieldSolver,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            /* first index is index of edge (x, y, z),
            second index is index of field component (ex, ey, ez or bx, by, bz) */
            const std::array<std::array<FunctionType, 3>, 3>& leftBFunc,
            const std::array<std::array<FunctionType, 3>, 3>& rightBFunc,
            const std::array<std::array<FunctionType, 3>, 3>& leftEFunc,
            const std::array<std::array<FunctionType, 3>, 3>& rightEFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1)) :
            FieldGenerator(fieldSolver, leftGenIndex, rightGenIndex,
                leftBFunc, rightBFunc, leftEFunc, rightEFunc,
                isLeftBorderEnabled, isRightBorderEnabled) {}

        void generateB() override {};
        void generateE() override {};
    };
}
