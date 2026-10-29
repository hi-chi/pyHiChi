#pragma once
#include <array>
#include <functional>

#include "Grid.h"
#include "Vectors.h"
#include "Enums.h"

namespace pfc
{
    template<class TGrid>
    class CurrentBoundaryCondition
    {
    public:

        CurrentBoundaryCondition(TGrid* grid):
            grid(grid),
            leftBorderIndex(grid->getNumExternalLeftCells()),
            rightBorderIndex(grid->getNumExternalLeftCells() + grid->numInternalCells)
        {}

        virtual ~CurrentBoundaryCondition() {}

        virtual void updateCurrentBoundaries() = 0;

        TGrid* grid = nullptr;
        Int3 leftBorderIndex, rightBorderIndex;
    };

    class PeriodicalCurrentBoundaryConditionFdtd : public CurrentBoundaryCondition<YeeGrid>
    {
    public:

        PeriodicalCurrentBoundaryConditionFdtd(YeeGrid* grid) :
            CurrentBoundaryCondition(grid)
        {}

        void updateCurrentBoundaries() override { 
            updateCurrent(pfc::CoordinateEnum::x);
            updateCurrent(pfc::CoordinateEnum::y);
            updateCurrent(pfc::CoordinateEnum::z);
        }

        void updateCurrent(CoordinateEnum axis);
    };

    inline void PeriodicalCurrentBoundaryConditionFdtd::updateCurrent(CoordinateEnum axis)
    {
        int dim0 = (int)axis;
        int dim1 = (dim0 + 1) % 3;
        int dim2 = (dim0 + 2) % 3;
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
                    indexR[dim0] = indexL[dim0] + grid->numInternalCells[dim0];
                    indexR[dim1] = j;
                    indexR[dim2] = k;

                    const int numExchangeCells = grid->getNumExternalLeftCells()[dim0];

                    for (int i = 0; i < numExchangeCells; i++) {
                        indexL[dim0] += i;
                        indexR[dim0] += i;

                        grid->Jx(indexR) += grid->Jx(indexL);
                        grid->Jy(indexR) += grid->Jy(indexL);
                        grid->Jz(indexR) += grid->Jz(indexL);

                        grid->Jx(indexL) = grid->Jx(indexR);
                        grid->Jy(indexL) = grid->Jy(indexR);
                        grid->Jz(indexL) = grid->Jz(indexR);
                    }

                    indexL[dim0] = grid->getNumExternalLeftCells()[dim0];
                    indexR[dim0] = indexL[dim0] + grid->numInternalCells[dim0];

                    for (int i = 0; i < numExchangeCells; i++) {
                        indexL[dim0] += i;
                        indexR[dim0] += i;

                        grid->Jx(indexL) += grid->Jx(indexR);
                        grid->Jy(indexL) += grid->Jy(indexR);
                        grid->Jz(indexL) += grid->Jz(indexR);

                        grid->Jx(indexR) = grid->Jx(indexL);
                        grid->Jy(indexR) = grid->Jy(indexL);
                        grid->Jz(indexR) = grid->Jz(indexL);
                    }
                }
    }
}
