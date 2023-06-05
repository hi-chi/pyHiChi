#pragma once
#include "FieldBoundaryCondition.h"

namespace pfc
{

    class PeriodicalBoundaryConditionFdtd : public FieldBoundaryCondition<GridTypes::YeeGridType>
    {
    public:
        PeriodicalBoundaryConditionFdtd(CoordinateEnum axis,
            FieldSolver<GridTypes::YeeGridType>* fieldSolver = 0) :
            FieldBoundaryCondition(axis, fieldSolver) {
        }

        PeriodicalBoundaryConditionFdtd(const PeriodicalBoundaryConditionFdtd& gen,
            FieldSolver<GridTypes::YeeGridType>* fieldSolver = 0) :
            FieldBoundaryCondition(gen, fieldSolver) {
        }

        void generateB() override;
        void generateE() override;

        FieldBoundaryCondition<GridTypes::YeeGridType>* createInstance(
            FieldSolver<GridTypes::YeeGridType>* fieldSolver) override {
            return new PeriodicalBoundaryConditionFdtd(*this, fieldSolver);
        }
    };

    inline void PeriodicalBoundaryConditionFdtd::generateB()
    {
        YeeGrid* grid = this->fieldSolver->grid;
        int dim0 = (int)axis;
        int dim1 = (dim0 + 1) % 3;
        int dim2 = (dim0 + 2) % 3;
        int begin1 = 0;
        int begin2 = 0;
        int end1 = grid->numCells[dim1];
        int end2 = grid->numCells[dim2];

        //OMP_FOR_COLLAPSE()
        for (int j = begin1; j < end1; j++)
            for (int k = begin2; k < end2; k++)
            {
                Int3 indexL, indexR;
                indexL[dim0] = 0;
                indexL[dim1] = j;
                indexL[dim2] = k;
                indexR[dim0] = indexL[dim0] + grid->numInternalCells[dim0];
                indexR[dim1] = j;
                indexR[dim2] = k;

                grid->Bx(indexL) = grid->Bx(indexR);
                grid->By(indexL) = grid->By(indexR);
                grid->Bz(indexL) = grid->Bz(indexR);
            }
    }

    inline void PeriodicalBoundaryConditionFdtd::generateE()
    {
        YeeGrid* grid = this->fieldSolver->grid;
        int dim0 = (int)axis;
        int dim1 = (dim0 + 1) % 3;
        int dim2 = (dim0 + 2) % 3;
        int begin1 = 0;
        int begin2 = 0;
        int end1 = grid->numCells[dim1];
        int end2 = grid->numCells[dim2];

        //OMP_FOR_COLLAPSE()
        for (int j = begin1; j < end1; j++)
            for (int k = begin2; k < end2; k++)
            {
                Int3 indexL, indexR;
                indexR[dim0] = grid->numCells[dim0] - 1;
                indexR[dim1] = j;
                indexR[dim2] = k;
                indexL[dim0] = indexR[dim0] - grid->numInternalCells[dim0];
                indexL[dim1] = j;
                indexL[dim2] = k;

                grid->Ex(indexR) = grid->Ex(indexL);
                grid->Ey(indexR) = grid->Ey(indexL);
                grid->Ez(indexR) = grid->Ez(indexL);
            }
    }


    class ReflectBoundaryConditionFdtd : public FieldBoundaryCondition<GridTypes::YeeGridType>
    {
    public:
        ReflectBoundaryConditionFdtd(CoordinateEnum axis,
            FieldSolver<GridTypes::YeeGridType>* fieldSolver = 0) :
            FieldBoundaryCondition(axis, fieldSolver) {
        }

        // copy constructor, other fieldSolver is possible
        ReflectBoundaryConditionFdtd(const FieldBoundaryCondition<GridTypes::YeeGridType>& gen,
            FieldSolver<GridTypes::YeeGridType>* fieldSolver = 0) :
            FieldBoundaryCondition(gen, fieldSolver) {
        }

        void generateB() override;
        void generateE() override;

        FieldBoundaryCondition<GridTypes::YeeGridType>* createInstance(
            FieldSolver<GridTypes::YeeGridType>* fieldSolver) override {
            return new ReflectBoundaryConditionFdtd(*this, fieldSolver);
        }
    };

    inline void ReflectBoundaryConditionFdtd::generateB()
    {
    }

    inline void ReflectBoundaryConditionFdtd::generateE()
    {
        YeeGrid* grid = this->fieldSolver->grid;
        int dim0 = (int)axis;
        int dim1 = (dim0 + 1) % 3;
        int dim2 = (dim0 + 2) % 3;
        int begin1 = this->fieldSolver->internalEAreaBegin[dim1];
        int begin2 = this->fieldSolver->internalEAreaBegin[dim2];
        int end1 = this->fieldSolver->internalEAreaEnd[dim1];
        int end2 = this->fieldSolver->internalEAreaEnd[dim2];

        //OMP_FOR_COLLAPSE()
        for (int j = begin1; j < end1; j++)
            for (int k = begin2; k < end2; k++)
            {
                // Adjust indexes for symmetry of generation coordinates
                Int3 indexL, indexR;
                indexL[dim0] = grid->getNumExternalLeftCells()[dim0] - 1;
                indexL[dim1] = j;
                indexL[dim2] = k;
                indexR[dim0] = indexL[dim0] + grid->numInternalCells[dim0];
                indexR[dim1] = j;
                indexR[dim2] = k;

                grid->Ex(indexL) = (FP)0.0;
                grid->Ey(indexL) = (FP)0.0;
                grid->Ez(indexL) = (FP)0.0;
                grid->Ex(indexR) = (FP)0.0;
                grid->Ey(indexR) = (FP)0.0;
                grid->Ez(indexR) = (FP)0.0;
            }
    }
}
