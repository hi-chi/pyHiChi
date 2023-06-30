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
            const Int3& domainIndexBegin, const Int3& domainIndexEnd,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
            FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGeneratorSpectral(TGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const std::array<FunctionType, 3>& xLeftBFunc,  // { bx, by, bz }
            const std::array<FunctionType, 3>& xRightBFunc,
            const std::array<FunctionType, 3>& yLeftBFunc,
            const std::array<FunctionType, 3>& yRightBFunc,
            const std::array<FunctionType, 3>& zLeftBFunc,
            const std::array<FunctionType, 3>& zRightBFunc,
            const std::array<FunctionType, 3>& xLeftEFunc,  // { ex, ey, ez }
            const std::array<FunctionType, 3>& xRightEFunc,
            const std::array<FunctionType, 3>& yLeftEFunc,
            const std::array<FunctionType, 3>& yRightEFunc,
            const std::array<FunctionType, 3>& zLeftEFunc,
            const std::array<FunctionType, 3>& zRightEFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGeneratorSpectral(TGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& bFunc,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& eFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGeneratorSpectral(TGrid* grid, FP dt, const Int3& domainIndexBegin,
            const Int3& domainIndexEnd, const FieldGeneratorSpectral& gen);

        // constructor for loading
        FieldGeneratorSpectral(TGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd);

        void generateB(FP time) {}
        void generateE(FP time) {}
    };

    template<class TGrid>
    FieldGeneratorSpectral<TGrid>::FieldGeneratorSpectral(TGrid* grid, FP dt,
        const Int3& domainIndexBegin, const Int3& domainIndexEnd,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
        FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd,
            leftGenIndex, rightGenIndex,
            bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
    {}

    template<class TGrid>
    FieldGeneratorSpectral<TGrid>::FieldGeneratorSpectral(TGrid* grid, FP dt,
        const Int3& domainIndexBegin, const Int3& domainIndexEnd,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<FunctionType, 3>& xLeftBFunc,  // { bx, by, bz }
        const std::array<FunctionType, 3>& xRightBFunc,
        const std::array<FunctionType, 3>& yLeftBFunc,
        const std::array<FunctionType, 3>& yRightBFunc,
        const std::array<FunctionType, 3>& zLeftBFunc,
        const std::array<FunctionType, 3>& zRightBFunc,
        const std::array<FunctionType, 3>& xLeftEFunc,  // { ex, ey, ez }
        const std::array<FunctionType, 3>& xRightEFunc,
        const std::array<FunctionType, 3>& yLeftEFunc,
        const std::array<FunctionType, 3>& yRightEFunc,
        const std::array<FunctionType, 3>& zLeftEFunc,
        const std::array<FunctionType, 3>& zRightEFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd,
            leftGenIndex, rightGenIndex,
            xLeftBFunc, xRightBFunc, yLeftBFunc, yRightBFunc, zLeftBFunc, zRightBFunc,
            xLeftEFunc, xRightEFunc, yLeftEFunc, yRightEFunc, zLeftEFunc, zRightEFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
    {}

    template<class TGrid>
    FieldGeneratorSpectral<TGrid>::FieldGeneratorSpectral(TGrid* grid, FP dt,
        const Int3& domainIndexBegin, const Int3& domainIndexEnd,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& bFunc,
        const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& eFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd,
            leftGenIndex, rightGenIndex,
            bFunc, eFunc, isLeftBorderEnabled, isRightBorderEnabled)
    {}

    template<class TGrid>
    FieldGeneratorSpectral<TGrid>::FieldGeneratorSpectral(TGrid* grid, FP dt, const Int3& domainIndexBegin,
        const Int3& domainIndexEnd, const FieldGeneratorSpectral& gen) :
        FieldGenerator<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd, gen)
    {}

    template<class TGrid>
    FieldGeneratorSpectral<TGrid>::FieldGeneratorSpectral(TGrid* grid, FP dt,
        const Int3& domainIndexBegin, const Int3& domainIndexEnd) :
        FieldGenerator<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd)
    {}
}
