#pragma once
#include "FieldBoundaryCondition.h"

namespace pfc
{

    class PeriodicalBoundaryConditionFdtd : public FieldBoundaryCondition<YeeGrid>
    {
    public:
        PeriodicalBoundaryConditionFdtd(YeeGrid* grid, CoordinateEnum axis) :
            FieldBoundaryCondition(grid, axis) {
        }

        void generateB(FP time) override;
        void generateE(FP time) override;

        FieldBoundaryCondition<YeeGrid>* createInstance(
            YeeGrid* grid, CoordinateEnum axis) override {
            return new PeriodicalBoundaryConditionFdtd(grid, axis);
        }
    };

    inline void PeriodicalBoundaryConditionFdtd::generateB(FP time)
    {
        int dim0 = (int)axis;
        int dim1 = (dim0 + 1) % 3;
        int dim2 = (dim0 + 2) % 3;
        // TODO: check border indices
        int begin1 = 0;
        int begin2 = 0;
        int end1 = this->grid->numCells[dim1];
        int end2 = this->grid->numCells[dim2];

        OMP_FOR_COLLAPSE()
        for (int j = begin1; j < end1; j++)
            for (int k = begin2; k < end2; k++)
            {
                Int3 indexL, indexR;
                indexL[dim0] = 0;
                indexL[dim1] = j;
                indexL[dim2] = k;
                indexR[dim0] = indexL[dim0] + this->grid->numInternalCells[dim0];
                indexR[dim1] = j;
                indexR[dim2] = k;

                this->grid->Bx(indexL) = this->grid->Bx(indexR);
                this->grid->By(indexL) = this->grid->By(indexR);
                this->grid->Bz(indexL) = this->grid->Bz(indexR);
            }
    }

    inline void PeriodicalBoundaryConditionFdtd::generateE(FP time)
    {
        int dim0 = (int)axis;
        int dim1 = (dim0 + 1) % 3;
        int dim2 = (dim0 + 2) % 3;
        // TODO: check border indices
        int begin1 = 0;
        int begin2 = 0;
        int end1 = this->grid->numCells[dim1];
        int end2 = this->grid->numCells[dim2];

        OMP_FOR_COLLAPSE()
        for (int j = begin1; j < end1; j++)
            for (int k = begin2; k < end2; k++)
            {
                Int3 indexL, indexR;
                indexR[dim0] = this->grid->numCells[dim0] - 1;
                indexR[dim1] = j;
                indexR[dim2] = k;
                indexL[dim0] = indexR[dim0] - this->grid->numInternalCells[dim0];
                indexL[dim1] = j;
                indexL[dim2] = k;

                this->grid->Ex(indexR) = this->grid->Ex(indexL);
                this->grid->Ey(indexR) = this->grid->Ey(indexL);
                this->grid->Ez(indexR) = this->grid->Ez(indexL);
            }
    }


    class ReflectBoundaryConditionFdtd : public FieldBoundaryCondition<YeeGrid>
    {
    public:
        ReflectBoundaryConditionFdtd(YeeGrid* grid, CoordinateEnum axis) :
            FieldBoundaryCondition(grid, axis) {
        }

        void generateB(FP time) override;
        void generateE(FP time) override;

        FieldBoundaryCondition<YeeGrid>* createInstance(
            YeeGrid* grid, CoordinateEnum axis) override {
            return new ReflectBoundaryConditionFdtd(grid, axis);
        }
    };

    inline void ReflectBoundaryConditionFdtd::generateB(FP time)
    {
    }

    inline void ReflectBoundaryConditionFdtd::generateE(FP time)
    {
        int dim0 = (int)axis;
        int dim1 = (dim0 + 1) % 3;
        int dim2 = (dim0 + 2) % 3;
        // TODO: check border indices
        int begin1 = 0;
        int begin2 = 0;
        int end1 = this->grid->numCells[dim1];
        int end2 = this->grid->numCells[dim2];

        OMP_FOR_COLLAPSE()
        for (int j = begin1; j < end1; j++)
            for (int k = begin2; k < end2; k++)
            {
                // Adjust indexes for symmetry of generation coordinates
                Int3 indexL, indexR;
                indexL[dim0] = this->grid->getNumExternalLeftCells()[dim0] - 1;
                indexL[dim1] = j;
                indexL[dim2] = k;
                indexR[dim0] = indexL[dim0] + this->grid->numInternalCells[dim0];
                indexR[dim1] = j;
                indexR[dim2] = k;

                this->grid->Ex(indexL) = (FP)0.0;
                this->grid->Ey(indexL) = (FP)0.0;
                this->grid->Ez(indexL) = (FP)0.0;
                this->grid->Ex(indexR) = (FP)0.0;
                this->grid->Ey(indexR) = (FP)0.0;
                this->grid->Ez(indexR) = (FP)0.0;
            }
    }
}
