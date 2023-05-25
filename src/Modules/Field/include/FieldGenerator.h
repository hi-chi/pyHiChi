#pragma once
#include <array>

#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"

namespace pfc
{
    template<GridTypes gridTypes>
    class FieldSolver;
    template<GridTypes gridTypes>
    class RealFieldSolver;
    template<GridTypes gridTypes>
    class SpectralFieldSolver;

    template<GridTypes gridTypes>
    class FieldGenerator
    {
    public:

        FieldGenerator(FieldSolver<gridTypes>* fieldSolver = 0);

        // copy constructor, other fieldSolver is possible
        FieldGenerator(const FieldGenerator& gen, FieldSolver<gridTypes>* fieldSolver = 0);

        virtual void generateB();
        virtual void generateE();

        FieldSolver<gridTypes>* fieldSolver;

        virtual FieldGenerator<gridTypes>* createInstance(FieldSolver<gridTypes>* fieldSolver) = 0;

    private:

        // major index is index of edge, minor index is index of component
        std::function<FP(FP, FP, FP, FP)> eLeft[3][3];
        std::function<FP(FP, FP, FP, FP)> eRight[3][3];
        std::function<FP(FP, FP, FP, FP)> bLeft[3][3];
        std::function<FP(FP, FP, FP, FP)> bRight[3][3];
        FP3 leftCoeff;
        FP3 rightCoeff;
    };

    typedef FieldGenerator<GridTypes::YeeGridType> FieldGeneratorYee;

    template<GridTypes gridTypes>
    inline FieldGenerator<gridTypes>::FieldGenerator(FieldSolver<gridTypes>* _fieldSolver)
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
                    /// ������, �� ����� ��������� ����� 1, � ���� �������� 1 ������ 2
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
}