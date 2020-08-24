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
        PSTD(PSTDGrid * grid);

        void updateFields();

        void updateHalfB();
        void updateE();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);

        void setTimeStep(FP dt);

    private:

        PmlSpectral<GridTypes::PSTDGridType>* getPml() {
            return (PmlSpectral<GridTypes::PSTDGridType>*)pml.get();
        }

    };

    inline PSTD::PSTD(PSTDGrid* grid) :
        SpectralFieldSolver<GridTypes::PSTDGridType>(grid)
    {
        updateDims();
        updateInternalDims();
    }

    inline void PSTD::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        pml.reset(new PmlPstd(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    inline void PSTD::setTimeStep(FP dt)
    {
        if (grid->setTimeStep(dt)) {
            complexGrid->setTimeStep(dt);
            if (pml.get()) pml.reset(new PmlPstd(this, pml->sizePML));
        }
    }

    inline void PSTD::updateFields()
    {
        doFourierTransform(fourier_transform::Direction::RtoC);

        if (pml.get()) getPml()->updateBSplit();
        updateHalfB();

        if (pml.get()) getPml()->updateESplit();
        updateE();

        if (pml.get()) getPml()->updateBSplit();
        updateHalfB();

        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml.get()) getPml()->doSecondStep();

        grid->globalTime += grid->dt;
        complexGrid->globalTime += grid->dt;
    }

    inline void PSTD::updateHalfB()
    {
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;
        double dt = grid->dt / 2;
#pragma omp parallel for collapse(2)
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
        double dt = grid->dt;
#pragma omp parallel for collapse(2)
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
