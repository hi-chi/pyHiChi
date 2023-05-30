#pragma once
#include <array>
#include <functional>

#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "Enums.h"

namespace pfc
{
    template<GridTypes gridTypes>
    class FieldSolver;

    namespace field_generator
    {
        inline FP defaultFieldFunction(FP x, FP y, FP z, FP t) {
            return (FP)0.0;
        }
    }

    template<GridTypes gridTypes>
    class FieldGenerator
    {
    public:

        using FunctionType = std::function<FP(FP, FP, FP, FP)>;  // f(x, y, z, t)

        FieldGenerator(FieldSolver<gridTypes>* fieldSolver,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGenerator(FieldSolver<gridTypes>* fieldSolver,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
            FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        FieldGenerator(FieldSolver<gridTypes>* fieldSolver,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            /* first index is index of edge (x, y, z),
            second index is index of field component (ex, ey, ez or bx, by, bz) */
            const std::array<std::array<FunctionType, 3>, 3>& leftBFunc,
            const std::array<std::array<FunctionType, 3>, 3>& rightBFunc,
            const std::array<std::array<FunctionType, 3>, 3>& leftEFunc,
            const std::array<std::array<FunctionType, 3>, 3>& rightEFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        virtual ~FieldGenerator() {}

        virtual void generateB() = 0;
        virtual void generateE() = 0;

        // sets one function for all borders
        void setFunction(FieldEnum field, CoordinateEnum fieldComponent, FunctionType func);
        // sets different functions for each border
        void setFunction(FieldEnum field, CoordinateEnum fieldComponent,
            CoordinateEnum edge, SideEnum side, FunctionType func);

        FieldSolver<gridTypes>* fieldSolver;

        // the first index is left/right
        // the second index is an index of edge
        // the third index is the index of field component
        std::array<std::array<std::array<FunctionType, 3>, 3>, 2> eFunc, bFunc;

        // sets enabled borders
        Int3 isLeftBorderEnabled, isRightBorderEnabled;
        Int3 leftGeneratorIndex, rightGeneratorIndex;
    };

    template<GridTypes gridTypes>
    inline FieldGenerator<gridTypes>::FieldGenerator(FieldSolver<gridTypes>* fieldSolver,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        fieldSolver(fieldSolver)
    {
        Int3 numExtLeftCells = fieldSolver->grid->getNumExternalLeftCells();

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

    template<GridTypes gridTypes>
    inline FieldGenerator<gridTypes>::FieldGenerator(FieldSolver<gridTypes>* fieldSolver,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
        FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<gridTypes>(fieldSolver, leftGenIndex, rightGenIndex,
            isLeftBorderEnabled, isRightBorderEnabled)
    {
        setFunction(FieldEnum::B, CoordinateEnum::x, bxFunc);
        setFunction(FieldEnum::B, CoordinateEnum::y, byFunc);
        setFunction(FieldEnum::B, CoordinateEnum::z, bzFunc);
        setFunction(FieldEnum::E, CoordinateEnum::x, exFunc);
        setFunction(FieldEnum::E, CoordinateEnum::y, eyFunc);
        setFunction(FieldEnum::E, CoordinateEnum::z, ezFunc);
    }

    template<GridTypes gridTypes>
    inline FieldGenerator<gridTypes>::FieldGenerator(FieldSolver<gridTypes>* fieldSolver,
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<std::array<FunctionType, 3>, 3>& leftBFunc,
        const std::array<std::array<FunctionType, 3>, 3>& rightBFunc,
        const std::array<std::array<FunctionType, 3>, 3>& leftEFunc,
        const std::array<std::array<FunctionType, 3>, 3>& rightEFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled) :
        FieldGenerator<gridTypes>(fieldSolver, leftGenIndex, rightGenIndex,
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

    template<GridTypes gridTypes>
    inline void FieldGenerator<gridTypes>::setFunction(
        FieldEnum field, CoordinateEnum fieldComponent,
        FieldGenerator<gridTypes>::FunctionType func)
    {
        for (int side = 0; side < 2; side++)
            for (int edge = 0; edge < 3; edge++)
                setFunction(field, fieldComponent,
                    (CoordinateEnum)edge, (SideEnum)side, func);
    }

    template<GridTypes gridTypes>
    inline void FieldGenerator<gridTypes>::setFunction(
        FieldEnum field, CoordinateEnum fieldComponent,
        CoordinateEnum edge, SideEnum side, FieldGenerator<gridTypes>::FunctionType func)
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


    template<GridTypes gridTypes>
    class SpectralFieldSolver;

    // temporary empty implementation of field generator for spectral field solvers
    template<GridTypes gridTypes>
    class SpectralFieldGenerator : public FieldGenerator<gridTypes>
    {
    public:

        SpectralFieldGenerator(FieldSolver<gridTypes>* fieldSolver,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
            FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1)) :
            FieldGenerator(fieldSolver, leftGenIndex, rightGenIndex,
                bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
                isLeftBorderEnabled, isRightBorderEnabled) {}

        SpectralFieldGenerator(FieldSolver<gridTypes>* fieldSolver,
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