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
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGenerator(TGrid* grid, FP dt,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
            FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGenerator(TGrid* grid, FP dt,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            /* first index is index of edge (x, y, z),
            second index is index of field component (ex, ey, ez or bx, by, bz) */
            const std::array<std::array<FunctionType, 3>, 3>& leftBFunc,
            const std::array<std::array<FunctionType, 3>, 3>& rightBFunc,
            const std::array<std::array<FunctionType, 3>, 3>& leftEFunc,
            const std::array<std::array<FunctionType, 3>, 3>& rightEFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGenerator(TGrid* grid, FP dt,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& bFunc,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& eFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGenerator(TGrid* grid, FP dt, const FieldGenerator& gen);

        /* implement the next methods in derived classes
        void generateB(FP time);
        void generateE(FP time);
        */

        // sets one function for all borders
        void setFunction(FieldEnum field, CoordinateEnum fieldComponent, FunctionType func);
        // sets different functions for each border
        void setFunction(FieldEnum field, CoordinateEnum fieldComponent,
            CoordinateEnum edge, SideEnum side, FunctionType func);

        TGrid* grid = nullptr;
        FP dt = 0.0;

        // first index is left/right
        // second index is index of edge
        // third index is index of field component
        std::array<std::array<std::array<FunctionType, 3>, 3>, 2> eFunc, bFunc;

        // sets enabled borders
        Int3 isLeftBorderEnabled, isRightBorderEnabled;
        Int3 leftGeneratorIndex, rightGeneratorIndex;
    };

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        grid(grid), dt(dt)
    {
        Int3 numExtLeftCells = this->grid->getNumExternalLeftCells();

        this->isLeftBorderEnabled = isLeftBorderEnabled;
        this->isRightBorderEnabled = isRightBorderEnabled;
        this->leftGeneratorIndex = leftGenIndex + numExtLeftCells;
        this->rightGeneratorIndex = rightGenIndex + numExtLeftCells;

        for (int side = 0; side < 2; side++)
            for (int edge = 0; edge < 3; edge++) {
                for (int comp = 0; comp < 3; comp++) {
                    eFunc[side][edge][comp] = field_generator::defaultFieldFunction;
                    bFunc[side][edge][comp] = field_generator::defaultFieldFunction;
                }
            }
    }

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
        FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<TGrid>(grid, dt, leftGenIndex, rightGenIndex,
            isLeftBorderEnabled, isRightBorderEnabled)
    {
        setFunction(FieldEnum::B, CoordinateEnum::x, bxFunc);
        setFunction(FieldEnum::B, CoordinateEnum::y, byFunc);
        setFunction(FieldEnum::B, CoordinateEnum::z, bzFunc);
        setFunction(FieldEnum::E, CoordinateEnum::x, exFunc);
        setFunction(FieldEnum::E, CoordinateEnum::y, eyFunc);
        setFunction(FieldEnum::E, CoordinateEnum::z, ezFunc);
    }

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<std::array<FunctionType, 3>, 3>& leftBFunc,
        const std::array<std::array<FunctionType, 3>, 3>& rightBFunc,
        const std::array<std::array<FunctionType, 3>, 3>& leftEFunc,
        const std::array<std::array<FunctionType, 3>, 3>& rightEFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<TGrid>(grid, dt, leftGenIndex, rightGenIndex,
            isLeftBorderEnabled, isRightBorderEnabled)
    {
        for (int edge = 0; edge < 3; edge++)
            for (int comp = 0; comp < 3; comp++) {
                setFunction(FieldEnum::B, (CoordinateEnum)comp,
                    (CoordinateEnum)edge, SideEnum::LEFT, leftBFunc[edge][comp]);
                setFunction(FieldEnum::B, (CoordinateEnum)comp,
                    (CoordinateEnum)edge, SideEnum::RIGHT, rightBFunc[edge][comp]);
                setFunction(FieldEnum::E, (CoordinateEnum)comp,
                    (CoordinateEnum)edge, SideEnum::LEFT, leftEFunc[edge][comp]);
                setFunction(FieldEnum::E, (CoordinateEnum)comp,
                    (CoordinateEnum)edge, SideEnum::RIGHT, rightEFunc[edge][comp]);
            }
    }

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& bFunc,
        const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& eFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<TGrid>(grid, dt, leftGenIndex, rightGenIndex,
            bFunc[0], bFunc[1], eFunc[0], eFunc[1],
            isLeftBorderEnabled, isRightBorderEnabled)
    {}

    template<class TGrid>
    inline FieldGenerator<TGrid>::FieldGenerator(TGrid* grid, FP dt,
        const FieldGenerator<TGrid>& gen) :
        FieldGenerator<TGrid>(grid, dt,
            gen.leftGeneratorIndex, gen.rightGeneratorIndex,
            gen.bFunc[0], gen.bFunc[1], gen.eFunc[0], gen.eFunc[1],
            gen.isLeftBorderEnabled, gen.isRightBorderEnabled)
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
}
