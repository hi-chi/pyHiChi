#pragma once
#include <array>

#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "Enums.h"

namespace pfc
{
    template<GridTypes gridTypes>
    class FieldSolver;

    namespace field_generator {
        inline FP defaultFieldFunction(FP x, FP y, FP z, FP t) {
            return (FP)0.0;
        }
    }

    template<GridTypes gridTypes>
    class FieldGenerator
    {
    public:

        typedef FP(*TFunc)(FP, FP, FP, FP);  // f(x, y, z, t)

        FieldGenerator(FieldSolver<gridTypes>* fieldSolver);

        void setFunction(FieldEnum field, CoordinateEnum edge, SideEnum side,
            CoordinateEnum fieldComponent, TFunc func);

        virtual void generateB() = 0;
        virtual void generateE() = 0;

        FieldSolver<gridTypes>* fieldSolver;

    protected:

        // the first index is left/right
        // the second index is an index of edge
        // the third index is the index of component
        std::array<std::array<std::array<TFunc, 3>, 3>, 2> eFunc, bFunc;

        // sets enabled borders
        std::array<std::array<bool, 3>, 2> isBorderEnabled;
        std::array<std::array<int, 3>, 2> generatorGridIndex;
    };

    template<GridTypes gridTypes>
    inline FieldGenerator<gridTypes>::FieldGenerator(FieldSolver<gridTypes>* fieldSolver) :
        fieldSolver(fieldSolver)
    {
        for (int size = 0; side < 2; side++)
            for (int edge = 0; edge < 3; edge++) {
                for (int comp = 0; comp < 3; comp++) {
                    eFunc[side][edge][comp] = field_generator::defaultFieldFunction;
                    bFunc[side][edge][comp] = field_generator::defaultFieldFunction;
                }
                isBorderEnabled[side][edge] = false;
                generatorGridIndex[side][edge] = -1;
            }
    }

    template<GridTypes gridTypes>
    inline void FieldGenerator<gridTypes>::setFunction(
        FieldEnum field, CoordinateEnum edge, SideEnum side,
        CoordinateEnum fieldComponent, FieldGenerator<gridTypes>::TFunc func)
    {
        switch (field) {
        case FieldEnum::B:
            bFunc[(int)side][(int)edge][(int)fieldComponent] = func;
            isBorderEnabled[(int)side][(int)edge] = true;
            break;
        case FieldEnum::E:
            eFunc[(int)side][(int)edge][(int)fieldComponent] = func;
            isBorderEnabled[(int)side][(int)edge] = true;
            break;
        case FieldEnum::J:
            std::cout << "WARNING: function for J is ignored in field generator" << std::endl;
            break;
        }
    }

    template<GridTypes gridTypes>
    inline void FieldGenerator<gridTypes>::generateB()
    {
        Grid<FP, gridTypes>* grid = fieldSolver->grid;
        const FP time = fieldSolver->globalTime +
            (grid->ifFieldsTimeStraggered ? fieldSolver->dt * 0.5 : 0.0);
        const FP cdt = constants::c * fieldSolver->dt;
        const FP3 norm_coeffs = FP3(cdt, cdt, cdt) / grid->steps;
        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
            int dim1 = (dim0 + 1) % 3;
            int dim2 = (dim0 + 2) % 3;
            int begin1 = fieldSolver->internalBAreaBegin[dim1];
            int begin2 = fieldSolver->internalBAreaBegin[dim2];
            int end1 = fieldSolver->internalBAreaEnd[dim1];
            int end2 = fieldSolver->internalBAreaEnd[dim2];
//OMP_FOR_COLLAPSE()
            for (int j = begin1; j < end1; j++)
                for (int k = begin2; k < end2; k++)
                {
                    // Adjust indexes for symmetry of generation coordinates
                    Int3 index;
                    index[dim0] = fieldSolver->internalBAreaBegin[dim0] + 1;
                    index[dim1] = j;
                    index[dim2] = k;
                    Int3 indexes[3] = { index, index, index };
                    indexes[dim1][dim0]++;
                    indexes[dim2][dim0]++;
                    FP3 bxCoords = grid->BxPosition(indexes[0].x, indexes[0].y,
                        indexes[0].z);
                    FP3 byCoords = grid->ByPosition(indexes[1].x, indexes[1].y,
                        indexes[1].z);
                    FP3 bzCoords = grid->BzPosition(indexes[2].x, indexes[2].y,
                        indexes[2].z);
                    FP coeff = leftCoeff[dim0] * norm_coeffs[dim0];
                    grid->Bx(indexes[0]) += coeff * bLeft[dim0][0](
                        bxCoords.x, bxCoords.y, bxCoords.z, time);
                    grid->By(indexes[1]) += coeff * bLeft[dim0][1](
                        byCoords.x, byCoords.y, byCoords.z, time);
                    grid->Bz(indexes[2]) += coeff * bLeft[dim0][2](
                        bzCoords.x, bzCoords.y, bzCoords.z, time);

                    index[dim0] = fieldSolver->internalBAreaEnd[dim0] - 2;
                    bxCoords = grid->BxPosition(index.x, index.y, index.z);
                    byCoords = grid->ByPosition(index.x, index.y, index.z);
                    bzCoords = grid->BzPosition(index.x, index.y, index.z);
                    coeff = rightCoeff[dim0] * norm_coeffs[dim0];
                    grid->Bx(index) += coeff * bRight[dim0][0](
                        bxCoords.x, bxCoords.y, bxCoords.z, time);
                    grid->By(index) += coeff * bRight[dim0][1](
                        byCoords.x, byCoords.y, byCoords.z, time);
                    grid->Bz(index) += coeff * bRight[dim0][2](
                        bzCoords.x, bzCoords.y, bzCoords.z, time);
                }
        }
    }

    template<GridTypes gridTypes>
    inline void FieldGenerator<gridTypes>::generateE()
    {
        Grid<FP, gridTypes> * grid = fieldSolver->grid;
        const FP time = fieldSolver->globalTime;
        const FP cdt = constants::c * fieldSolver->dt;
        const FP3 norm_coeffs = FP3(cdt, cdt, cdt) / grid->steps;
        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
            int dim1 = (dim0 + 1) % 3;
            int dim2 = (dim0 + 2) % 3;
            int begin1 = fieldSolver->internalEAreaBegin[dim1];
            int begin2 = fieldSolver->internalEAreaBegin[dim2];
            int end1 = fieldSolver->internalEAreaEnd[dim1];
            int end2 = fieldSolver->internalEAreaEnd[dim2];
//OMP_FOR_COLLAPSE()
            for (int j = begin1; j < end1; j++)
                for (int k = begin2; k < end2; k++)
                {
                    // Adjust indexes for symmetry of generation coordinates, use B indexes for consistency
                    /// Видимо, не нужно добавлять здесь 1, а ниже вычитать 1 вместо 2
                    Int3 index;
                    index[dim0] = fieldSolver->internalBAreaBegin[dim0] + 1;
                    index[dim1] = j;
                    index[dim2] = k;
                    Int3 indexes[3] = { index, index, index };
                    indexes[dim0][dim0]++;
                    FP3 exCoords = grid->ExPosition(indexes[0].x, indexes[0].y,
                        indexes[0].z);
                    FP3 eyCoords = grid->EyPosition(indexes[1].x, indexes[1].y,
                        indexes[1].z);
                    FP3 ezCoords = grid->EzPosition(indexes[2].x, indexes[2].y,
                        indexes[2].z);
                    FP coeff = leftCoeff[dim0] * norm_coeffs[dim0];
                    grid->Ex(indexes[0]) += coeff * eLeft[dim0][0](
                        exCoords.x, exCoords.y, exCoords.z, time);
                    grid->Ey(indexes[1]) += coeff * eLeft[dim0][1](
                        eyCoords.x, eyCoords.y, eyCoords.z, time);
                    grid->Ez(indexes[2]) += coeff * eLeft[dim0][2](
                        ezCoords.x, ezCoords.y, ezCoords.z, time);

                    index[dim0] = fieldSolver->internalBAreaEnd[dim0] - 2;
                    exCoords = grid->ExPosition(index.x, index.y, index.z);
                    eyCoords = grid->EyPosition(index.x, index.y, index.z);
                    ezCoords = grid->EzPosition(index.x, index.y, index.z);
                    coeff = rightCoeff[dim0] * norm_coeffs[dim0];
                    grid->Ex(index) += coeff * eRight[dim0][0](
                        exCoords.x, exCoords.y, exCoords.z, time);
                    grid->Ey(index) += coeff * eRight[dim0][1](
                        eyCoords.x, eyCoords.y, eyCoords.z, time);
                    grid->Ez(index) += coeff * eRight[dim0][2](
                        ezCoords.x, ezCoords.y, ezCoords.z, time);
                }
        }
    }

    typedef FieldGenerator<GridTypes::YeeGridType> FieldGeneratorYee;
}