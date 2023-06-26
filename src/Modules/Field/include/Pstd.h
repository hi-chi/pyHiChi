#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "PmlPstd.h"
#include "FieldBoundaryConditionSpectral.h"
#include "FieldGeneratorSpectral.h"

namespace pfc {
    class PSTD : public SpectralFieldSolver<PSTDGrid, PmlPstd, FieldGeneratorSpectral<PSTDGrid>>
    {
    public:

        using GridType = PSTDGrid;
        using PmlType = PmlPstd;
        using FieldGeneratorType = FieldGeneratorSpectral<PSTDGrid>;
        using PeriodicalBoundaryConditionType = PeriodicalBoundaryConditionSpectral<PSTDGrid>;

        PSTD(GridType* grid);  // use when load
        PSTD(GridType* grid, FP dt);

        void updateFields();

        void updateHalfB();
        void updateE();

        void setTimeStep(FP dt);

        static FP getCourantConditionTimeStep(const FP3& gridSteps) {
            FP tmp = sqrt(1.0 / (gridSteps.x * gridSteps.x) +
                1.0 / (gridSteps.y * gridSteps.y) +
                1.0 / (gridSteps.z * gridSteps.z));
            return 2.0 / (constants::pi * constants::c * tmp);
        }
        FP getCourantConditionTimeStep() const {
            return getCourantConditionTimeStep(grid->steps);
        }

        bool isCourantConditionSatisfied(FP dt) const {
            return dt < getCourantConditionTimeStep();
        }

        static bool isTimeStaggered() {
            return true;
        }
    };

    inline PSTD::PSTD(GridType* grid) :
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>(grid)
    {}

    inline PSTD::PSTD(GridType* grid, FP dt) :
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>(grid, dt)
    {
        if (!isCourantConditionSatisfied(dt)) {
            std::cout
                << "WARNING: PSTD Courant condition is not satisfied. Another time step was setted up"
                << std::endl;
            this->dt = getCourantConditionTimeStep() * 0.5;
        }
    }

    inline void PSTD::setTimeStep(FP dt)
    {
        if (isCourantConditionSatisfied(dt)) {
            this->dt = dt;
            resetPML();
            resetFieldGenerator();
        }
        else {
            std::cout
                << "WARNING: PSTD Courant condition is not satisfied. Time step was not changed"
                << std::endl;
        }
    }

    inline void PSTD::updateFields()
    {
        // TODO: consider boundary conditions and generator

        doFourierTransform(fourier_transform::Direction::RtoC);

        if (pml) pml->updateBSplit();
        updateHalfB();
        // if (generator) generator->generateB(globalTime);  // send current E time
        // applyBoundaryConditionsB(globalTime + dt * 0.5);

        if (pml) pml->updateESplit();
        updateE();
        // if (generator) generator->generateE(globalTime + dt * 0.5);  // send current B time
        // applyBoundaryConditionsE(globalTime + dt);

        updateHalfB();
        // applyBoundaryConditionsB(globalTime + dt);

        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml) pml->updateB();
        if (pml) pml->updateE();

        globalTime += dt;
    }

    inline void PSTD::updateHalfB()
    {
        const Int3 begin = this->complexDomainIndexBegin;
        const Int3 end = this->complexDomainIndexEnd;

        FP dt = 0.5 * this->dt;

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 crossKE = cross((ComplexFP3)getWaveVector(Int3(i, j, k)), E);
                    complexFP coeff = -complexFP::i() * constants::c * dt;

                    complexGrid->Bx(i, j, k) += coeff * crossKE.x;
                    complexGrid->By(i, j, k) += coeff * crossKE.y;
                    complexGrid->Bz(i, j, k) += coeff * crossKE.z;
                }
    }

    inline void PSTD::updateE()
    {
        const Int3 begin = this->complexDomainIndexBegin;
        const Int3 end = this->complexDomainIndexEnd;

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
                    ComplexFP3 B(complexGrid->Bx(i, j, k), complexGrid->By(i, j, k), complexGrid->Bz(i, j, k));
                    ComplexFP3 J(complexGrid->Jx(i, j, k), complexGrid->Jy(i, j, k), complexGrid->Jz(i, j, k));
                    ComplexFP3 crossKB = cross((ComplexFP3)getWaveVector(Int3(i, j, k)), B);
                    complexFP coeff = complexFP::i() * constants::c * dt;

                    complexGrid->Ex(i, j, k) += coeff * crossKB.x - 4 * constants::pi * dt * J.x;
                    complexGrid->Ey(i, j, k) += coeff * crossKB.y - 4 * constants::pi * dt * J.y;
                    complexGrid->Ez(i, j, k) += coeff * crossKB.z - 4 * constants::pi * dt * J.z;
                }
    }

}
