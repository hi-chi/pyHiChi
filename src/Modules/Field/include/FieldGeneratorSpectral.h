#pragma once
#include "FieldGenerator.h"

namespace pfc
{
    // temporary empty implementation of field generator for spectral field solvers
    template<class TGrid>
    class FieldGeneratorSpectral : public FieldGenerator<TGrid>
    {
    public:
    
        using FunctionType = typename FieldGenerator<TGrid>::FunctionType;

        FieldGeneratorSpectral(TGrid* grid, FP dt,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
            FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1)) :
            FieldGenerator<TGrid>(grid, dt, leftGenIndex, rightGenIndex,
                bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
                isLeftBorderEnabled, isRightBorderEnabled) {}

        FieldGeneratorSpectral(TGrid* grid, FP dt,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            /* first index is index of edge (x, y, z),
            second index is index of field component (ex, ey, ez or bx, by, bz) */
            const std::array<std::array<FunctionType, 3>, 3>& leftBFunc,
            const std::array<std::array<FunctionType, 3>, 3>& rightBFunc,
            const std::array<std::array<FunctionType, 3>, 3>& leftEFunc,
            const std::array<std::array<FunctionType, 3>, 3>& rightEFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1)) :
            FieldGenerator<TGrid>(grid, dt, leftGenIndex, rightGenIndex,
                leftBFunc, rightBFunc, leftEFunc, rightEFunc,
                isLeftBorderEnabled, isRightBorderEnabled) {}

        FieldGeneratorSpectral(TGrid* grid, FP dt,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& bFunc,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& eFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1)) :
            FieldGenerator<TGrid>(grid, dt, leftGenIndex, rightGenIndex,
                bFunc, eFunc, isLeftBorderEnabled, isRightBorderEnabled) {}

        FieldGeneratorSpectral(TGrid* grid, FP dt, const FieldGeneratorSpectral& gen) :
            FieldGenerator<TGrid>(grid, dt, gen) {}

        void generateB(FP time) {}
        void generateE(FP time) {}
    };

}
