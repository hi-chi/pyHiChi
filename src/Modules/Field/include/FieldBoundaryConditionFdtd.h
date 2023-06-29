#pragma once
#include "FieldBoundaryCondition.h"

namespace pfc
{

    class PeriodicalBoundaryConditionFdtd : public FieldBoundaryCondition<YeeGrid>
    {
    public:

        PeriodicalBoundaryConditionFdtd(YeeGrid* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) :
            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex, axis)
        {}

        // constructor for loading
        PeriodicalBoundaryConditionFdtd(YeeGrid* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex) :
            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex)
        {}

        void generateB(FP time) override { generateField(time, grid->Bx, grid->By, grid->Bz); }
        void generateE(FP time) override { generateField(time, grid->Ex, grid->Ey, grid->Ez); }

        void generateField(FP time, ScalarField<FP>& fx, ScalarField<FP>& fy, ScalarField<FP>& fz);

        FieldBoundaryCondition<YeeGrid>* createInstance(
            YeeGrid* grid, Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) override {
            return new PeriodicalBoundaryConditionFdtd(grid, leftBorderIndex, rightBorderIndex, axis);
        }
    };

    inline void PeriodicalBoundaryConditionFdtd::generateField(FP time,
        ScalarField<FP>& fx, ScalarField<FP>& fy, ScalarField<FP>& fz)
    {
        int dim0 = (int)axis;
        int dim1 = (dim0 + 1) % 3;
        int dim2 = (dim0 + 2) % 3;
        int begin1 = this->leftBorderIndex[dim1];
        int begin2 = this->leftBorderIndex[dim2];
        int end1 = this->rightBorderIndex[dim1];
        int end2 = this->rightBorderIndex[dim2];

        OMP_FOR_COLLAPSE()
        for (int j = begin1; j < end1; j++)
            for (int k = begin2; k < end2; k++)
            {
                Int3 indexL, indexR;
                indexL[dim1] = indexR[dim1] = j;
                indexL[dim2] = indexR[dim2] = k;

                indexL[dim0] = this->leftBorderIndex[dim0] - 1;
                indexR[dim0] = this->rightBorderIndex[dim0] - 1;

                fx(indexL) = fx(indexR);
                fy(indexL) = fy(indexR);
                fz(indexL) = fz(indexR);

                indexL[dim0]++;
                indexR[dim0]++;

                fx(indexR) = fx(indexL);
                fy(indexR) = fy(indexL);
                fz(indexR) = fz(indexL);
            }
    }


    class ReflectBoundaryConditionFdtd : public FieldBoundaryCondition<YeeGrid>
    {
    public:

        ReflectBoundaryConditionFdtd(YeeGrid* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) :
            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex, axis)
        {}

        // constructor for loading
        ReflectBoundaryConditionFdtd(YeeGrid* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex) :
            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex)
        {}

        void generateB(FP time) override {}
        void generateE(FP time) override;

        FieldBoundaryCondition<YeeGrid>* createInstance(
            YeeGrid* grid, Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) override {
            return new ReflectBoundaryConditionFdtd(grid, leftBorderIndex, rightBorderIndex, axis);
        }
    };

    inline void ReflectBoundaryConditionFdtd::generateE(FP time)
    {
        int dim0 = (int)axis;
        int dim1 = (dim0 + 1) % 3;
        int dim2 = (dim0 + 2) % 3;
        int begin1 = this->leftBorderIndex[dim1];
        int begin2 = this->leftBorderIndex[dim2];
        int end1 = this->rightBorderIndex[dim1];
        int end2 = this->rightBorderIndex[dim2];

        OMP_FOR_COLLAPSE()
        for (int j = begin1; j < end1; j++)
            for (int k = begin2; k < end2; k++)
            {
                Int3 indexL, indexR;
                indexL[dim1] = indexR[dim1] = j;
                indexL[dim2] = indexR[dim2] = k;

                indexL[dim0] = this->leftBorderIndex[dim0] - 1;
                indexR[dim0] = this->rightBorderIndex[dim0] - 1;

                this->grid->Ex(indexL) = (FP)0.0;
                this->grid->Ey(indexL) = (FP)0.0;
                this->grid->Ez(indexL) = (FP)0.0;
                this->grid->Ex(indexR) = (FP)0.0;
                this->grid->Ey(indexR) = (FP)0.0;
                this->grid->Ez(indexR) = (FP)0.0;
            }
    }
}
