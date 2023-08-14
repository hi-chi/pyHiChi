#pragma once
#include <array>
#include <functional>

#include "Grid.h"
#include "Vectors.h"
#include "Enums.h"

namespace pfc
{
    namespace field_generator
    {
        inline FP defaultFieldFunction(FP x, FP y, FP z, FP t) {
            return (FP)0.0;
        }
    }

    template<class TGrid>
    class FieldGenerator
    {
    public:

        using FunctionType = std::function<FP(FP, FP, FP, FP)>;  // f(x, y, z, t)

        FieldGenerator(TGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
            FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGenerator(TGrid* grid, FP dt,
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

        FieldGenerator(TGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& bFunc,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& eFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGenerator(TGrid* grid, FP dt, Int3 domainIndexBegin,
            Int3 domainIndexEnd, const FieldGenerator& gen);

        // constructor for loading
        explicit FieldGenerator(TGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd);

        /* implement the next methods in derived classes
        void generateB(FP time);
        void generateE(FP time);
        */

        // sets one function for all borders
        void setFunction(FieldEnum field, CoordinateEnum fieldComponent, FunctionType func);
        // sets different functions for each border
        void setFunction(FieldEnum field, CoordinateEnum fieldComponent,
            CoordinateEnum edge, SideEnum side, FunctionType func);

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        TGrid* grid = nullptr;
        FP dt = 0.0;
        Int3 domainIndexBegin, domainIndexEnd;

        // generator position after domainIndexBegin
        Int3 leftGeneratorIndex, rightGeneratorIndex;
        // sets enabled borders
        Int3 isLeftBorderEnabled, isRightBorderEnabled;

        // first index is left/right
        // second index is index of edge
        // third index is index of field component
        std::array<std::array<std::array<FunctionType, 3>, 3>, 2> eFunc, bFunc;
    };

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        const Int3& domainIndexBegin, const Int3& domainIndexEnd,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
        FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd,
            leftGenIndex, rightGenIndex,
            { bxFunc, byFunc, bzFunc }, { bxFunc, byFunc, bzFunc },
            { bxFunc, byFunc, bzFunc }, { bxFunc, byFunc, bzFunc },
            { bxFunc, byFunc, bzFunc }, { bxFunc, byFunc, bzFunc },
            { exFunc, eyFunc, ezFunc }, { exFunc, eyFunc, ezFunc },
            { exFunc, eyFunc, ezFunc }, { exFunc, eyFunc, ezFunc },
            { exFunc, eyFunc, ezFunc }, { exFunc, eyFunc, ezFunc },
            isLeftBorderEnabled, isRightBorderEnabled)
    {}

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        const Int3& domainIndexBegin, const Int3& domainIndexEnd,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<FunctionType, 3>& xLeftBFunc,
        const std::array<FunctionType, 3>& xRightBFunc,
        const std::array<FunctionType, 3>& yLeftBFunc,
        const std::array<FunctionType, 3>& yRightBFunc,
        const std::array<FunctionType, 3>& zLeftBFunc,
        const std::array<FunctionType, 3>& zRightBFunc,
        const std::array<FunctionType, 3>& xLeftEFunc,
        const std::array<FunctionType, 3>& xRightEFunc,
        const std::array<FunctionType, 3>& yLeftEFunc,
        const std::array<FunctionType, 3>& yRightEFunc,
        const std::array<FunctionType, 3>& zLeftEFunc,
        const std::array<FunctionType, 3>& zRightEFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd,
            leftGenIndex, rightGenIndex,
            { std::array<std::array<FunctionType, 3>, 3>{ xLeftBFunc, yLeftBFunc, zLeftBFunc },
              std::array<std::array<FunctionType, 3>, 3>{ xRightBFunc, yRightBFunc, zRightBFunc } },
            { std::array<std::array<FunctionType, 3>, 3>{ xLeftEFunc, yLeftEFunc, zLeftEFunc },
              std::array<std::array<FunctionType, 3>, 3>{ xRightEFunc, yRightEFunc, zRightEFunc } },
            isLeftBorderEnabled, isRightBorderEnabled)
    {}

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        const Int3& domainIndexBegin, const Int3& domainIndexEnd,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& bFunc,
        const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& eFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        grid(grid), dt(dt),
        domainIndexBegin(domainIndexBegin), domainIndexEnd(domainIndexEnd),
        leftGeneratorIndex(domainIndexBegin + leftGenIndex),
        rightGeneratorIndex(domainIndexBegin + rightGenIndex),
        bFunc(bFunc), eFunc(eFunc),
        isLeftBorderEnabled(isLeftBorderEnabled), isRightBorderEnabled(isRightBorderEnabled)
    {}

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        Int3 domainIndexBegin, Int3 domainIndexEnd, const FieldGenerator<TGrid>& gen) :
        FieldGenerator<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd,
            gen.leftGeneratorIndex, gen.rightGeneratorIndex,
            gen.bFunc, gen.eFunc,
            gen.isLeftBorderEnabled, gen.isRightBorderEnabled)
    {}

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        const Int3& domainIndexBegin, const Int3& domainIndexEnd) :
        grid(grid), dt(dt),
        domainIndexBegin(domainIndexBegin), domainIndexEnd(domainIndexEnd)
    {}

    template<class TGrid>
    inline void FieldGenerator<TGrid>::setFunction(
        FieldEnum field, CoordinateEnum fieldComponent,
        FieldGenerator<TGrid>::FunctionType func)
    {
        for (int side = 0; side < 2; side++)
            for (int edge = 0; edge < 3; edge++)
                setFunction(field, fieldComponent,
                    (CoordinateEnum)edge, (SideEnum)side, func);
    }

    template<class TGrid>
    inline void FieldGenerator<TGrid>::setFunction(
        FieldEnum field, CoordinateEnum fieldComponent,
        CoordinateEnum edge, SideEnum side, FieldGenerator<TGrid>::FunctionType func)
    {
        switch (field) {
        case FieldEnum::B:
            bFunc[(int)side][(int)edge][(int)fieldComponent] = func;
            break;
        case FieldEnum::E:
            eFunc[(int)side][(int)edge][(int)fieldComponent] = func;
            break;
        default:
            break;
        }
    }

    template<class TGrid>
    inline void FieldGenerator<TGrid>::save(std::ostream& ostr)
    {
        ostr.write((char*)&leftGeneratorIndex, sizeof(leftGeneratorIndex));
        ostr.write((char*)&rightGeneratorIndex, sizeof(rightGeneratorIndex));
        ostr.write((char*)&isLeftBorderEnabled, sizeof(isLeftBorderEnabled));
        ostr.write((char*)&isRightBorderEnabled, sizeof(isRightBorderEnabled));
        // TODO: save functions
    }

    template<class TGrid>
    inline void FieldGenerator<TGrid>::load(std::istream& istr)
    {
        istr.read((char*)&leftGeneratorIndex, sizeof(leftGeneratorIndex));
        istr.read((char*)&rightGeneratorIndex, sizeof(rightGeneratorIndex));
        istr.read((char*)&isLeftBorderEnabled, sizeof(isLeftBorderEnabled));
        istr.read((char*)&isRightBorderEnabled, sizeof(isRightBorderEnabled));
        // TODO: load functions
    }
}
