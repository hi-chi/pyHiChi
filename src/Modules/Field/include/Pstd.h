#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "PmlPstd.h"

namespace pfc {
    class PSTD : public SpectralFieldSolver<PSTDGridType>
    {
    public:

        using GridType = PSTDGrid;
        using PmlType = PmlPstd;
        using PeriodicalFieldGeneratorType = PeriodicalFieldGeneratorPstd;

        PSTD(GridType* grid, double dt);

        void updateFields();

        void updateHalfB();
        void updateE();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void setFieldGenerator(PeriodicalFieldGeneratorType* _generator);

        void setTimeStep(FP dt);

        FP getCourantCondition() const {
            double tmp = sqrt(1.0 / (grid->steps.x*grid->steps.x) +
                1.0 / (grid->steps.y*grid->steps.y) +
                1.0 / (grid->steps.z*grid->steps.z));
            return 2.0 / (constants::pi * constants::c * tmp);
        }

        bool ifCourantConditionSatisfied(FP dt) const {
            return dt < getCourantCondition();
        }

    private:

        PmlSpectral<GridTypes::PSTDGridType>* getPml() {
            return (PmlSpectral<GridTypes::PSTDGridType>*)pml.get();
        }

    };

    inline PSTD::PSTD(PSTD::GridType* grid, FP dt) :
        SpectralFieldSolver(grid, dt, 0.0, 0.5 * dt, 0.5 * dt)
    {
        if (!ifCourantConditionSatisfied(dt)) {
            std::cout
                << "WARNING: PSTD Courant condition is not satisfied. Another time step was setted up"
                << std::endl;
            this->dt = getCourantCondition() * 0.5;
        }
        updateDims();
        updateInternalDims();
    }

    inline void PSTD::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        pml.reset(new PmlType(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    inline void PSTD::setFieldGenerator(PSTD::PeriodicalFieldGeneratorType* _generator)
    {
        generator.reset(_generator->createInstance(this));
    }

    inline void PSTD::setTimeStep(FP dt)
    {
        if (ifCourantConditionSatisfied(dt)) {
            this->dt = dt;
            this->timeShiftB = 0.5 * dt;
            this->timeShiftJ = 0.5 * dt;
            if (pml) pml.reset(new PmlType(this, pml->sizePml));
            if (generator) generator.reset(generator->createInstance(this));
        }
        else {
            std::cout
                << "WARNING: PSTD Courant condition is not satisfied. Time step was not changed"
                << std::endl;
        }
    }

    inline void PSTD::updateFields()
    {
        doFourierTransform(fourier_transform::Direction::RtoC);

        if (pml) getPml()->updateBSplit();
        updateHalfB();

        if (pml) getPml()->updateESplit();
        updateE();

        if (pml) getPml()->updateBSplit();
        updateHalfB();

        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml) getPml()->updateB();
        if (pml) getPml()->updateE();

        globalTime += dt;
    }

    inline void PSTD::updateHalfB()
    {
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;
        double dt = 0.5 * this->dt;
        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
            {
//#pragma omp simd
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
    }

    inline void PSTD::updateE()
    {
        const Int3 begin = updateComplexEAreaBegin;
        const Int3 end = updateComplexEAreaEnd;
        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
            {
//#pragma omp simd
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

}
