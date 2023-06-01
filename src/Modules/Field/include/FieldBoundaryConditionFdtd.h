#pragma once
#include "FieldBoundaryCondition.h"

namespace pfc
{

    class PeriodicalBoundaryConditionFdtd : public FieldBoundaryCondition<GridTypes::YeeGridType>
    {
    public:
        PeriodicalBoundaryConditionFdtd(FieldSolver<GridTypes::YeeGridType>* fieldSolver = 0,
            bool isXAxisEnabled = true, bool isYAxisEnabled = true, bool isZAxisEnabled = true) :
            FieldBoundaryCondition(fieldSolver, isXAxisEnabled, isYAxisEnabled, isZAxisEnabled) {
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
        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
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
                    indexL[dim0] = this->fieldSolver->internalEAreaBegin[dim0];
                    indexL[dim1] = j;
                    indexL[dim2] = k;
                    indexR[dim0] = this->fieldSolver->internalEAreaEnd[dim0] - 3;
                    indexR[dim1] = j;
                    indexR[dim2] = k;

                    grid->Bx(indexL) = grid->Bx(indexR);
                    grid->By(indexL) = grid->By(indexR);
                    grid->Bz(indexL) = grid->Bz(indexR);

                    indexL[dim0]++;
                    indexR[dim0]++;

                    grid->Bx(indexR) = grid->Bx(indexL);
                    grid->By(indexR) = grid->By(indexL);
                    grid->Bz(indexR) = grid->Bz(indexL);
                }
        }
    }

    inline void PeriodicalBoundaryConditionFdtd::generateE()
    {
        YeeGrid* grid = this->fieldSolver->grid;
        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
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
                    // Adjust indexes for symmetry of generation coordinates
                    Int3 indexL, indexR;
                    indexL[dim0] = this->fieldSolver->internalEAreaBegin[dim0] + 1;
                    indexL[dim1] = j;
                    indexL[dim2] = k;
                    indexR[dim0] = this->fieldSolver->internalEAreaEnd[dim0] - 2;
                    indexR[dim1] = j;
                    indexR[dim2] = k;

                    grid->Ex(indexL) = grid->Ex(indexR);
                    grid->Ey(indexL) = grid->Ey(indexR);
                    grid->Ez(indexL) = grid->Ez(indexR);

                    indexL[dim0]++;
                    indexR[dim0]++;

                    grid->Ex(indexR) = grid->Ex(indexL);
                    grid->Ey(indexR) = grid->Ey(indexL);
                    grid->Ez(indexR) = grid->Ez(indexL);
                }
        }
    }


    class ReflectBoundaryConditionFdtd : public FieldBoundaryCondition<GridTypes::YeeGridType>
    {
    public:
        ReflectBoundaryConditionFdtd(FieldSolver<GridTypes::YeeGridType>* fieldSolver = 0,
            bool isXAxisEnabled = true, bool isYAxisEnabled = true, bool isZAxisEnabled = true) :
            FieldBoundaryCondition(
                fieldSolver, isXAxisEnabled, isYAxisEnabled, isZAxisEnabled) {
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
        YeeGrid* grid = this->fieldSolver->grid;
        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
            int dim1 = (dim0 + 1) % 3;
            int dim2 = (dim0 + 2) % 3;
            int begin1 = this->fieldSolver->internalBAreaBegin[dim1];
            int begin2 = this->fieldSolver->internalBAreaBegin[dim2];
            int end1 = this->fieldSolver->internalBAreaEnd[dim1];
            int end2 = this->fieldSolver->internalBAreaEnd[dim2];
            //OMP_FOR_COLLAPSE()
            for (int j = begin1; j < end1; j++)
                for (int k = begin2; k < end2; k++)
                {
                    // Adjust indexes for symmetry of generation coordinates
                    Int3 indexL, indexR;
                    indexL[dim0] = this->fieldSolver->internalBAreaBegin[dim0];
                    indexL[dim1] = j;
                    indexL[dim2] = k;
                    indexR[dim0] = this->fieldSolver->internalBAreaEnd[dim0] - 1;
                    indexR[dim1] = j;
                    indexR[dim2] = k;

                    grid->Bx(indexR) = (FP)0.0;
                    grid->By(indexR) = (FP)0.0;
                    grid->Bz(indexR) = (FP)0.0;
                    grid->Bx(indexL) = (FP)0.0;
                    grid->By(indexL) = (FP)0.0;
                    grid->Bz(indexL) = (FP)0.0;
                }
        }
    }

    inline void ReflectBoundaryConditionFdtd::generateE()
    {
        YeeGrid* grid = this->fieldSolver->grid;
        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
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
                    indexL[dim0] = this->fieldSolver->internalEAreaBegin[dim0];
                    indexL[dim1] = j;
                    indexL[dim2] = k;
                    indexR[dim0] = this->fieldSolver->internalEAreaEnd[dim0] - 1;
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
}
