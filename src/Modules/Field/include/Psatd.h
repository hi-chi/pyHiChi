#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "PmlPsatd.h"
//#include <chrono>
#include <omp.h>

namespace pfc {

    class PSATDTimeStraggered : public SpectralFieldSolver<PSATDTimeStraggeredGridType>
    {

    public:
        PSATDTimeStraggered(PSATDTimeStraggeredGrid * grid);

        void updateFields();

        void updateHalfB();
        void updateE();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);

        void setTimeStep(FP dt);

        void convertFieldsPoissonEquation();

        ScalarField<complexFP> tmpJx, tmpJy, tmpJz;

    protected:

        PmlSpectral<GridTypes::PSATDTimeStraggeredGridType>* getPml() {
            return (PmlSpectral<GridTypes::PSATDTimeStraggeredGridType>*)pml.get();
        }

        void saveJ();
        void assignJ(ScalarField<complexFP>& J, ScalarField<complexFP>& tmpJ);
    };

    inline PSATDTimeStraggered::PSATDTimeStraggered(PSATDTimeStraggeredGrid* _grid) :
        SpectralFieldSolver<GridTypes::PSATDTimeStraggeredGridType>(_grid),
        tmpJx(FourierTransform::getSizeOfComplex(_grid->Jx.getSize())),
        tmpJy(FourierTransform::getSizeOfComplex(_grid->Jy.getSize())),
        tmpJz(FourierTransform::getSizeOfComplex(_grid->Jz.getSize()))
    {
        updateDims();
        updateInternalDims();
    }

    inline void PSATDTimeStraggered::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        pml.reset(new PmlPsatdTimeStraggered(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    inline void PSATDTimeStraggered::setTimeStep(FP dt)
    {
        if (grid->setTimeStep(dt)) {
            complexGrid->setTimeStep(dt);
            if (pml.get()) pml.reset(new PmlPsatdTimeStraggered(this, pml->sizePML));
        }
    }

    inline void PSATDTimeStraggered::assignJ(ScalarField<complexFP>& J, ScalarField<complexFP>& tmpJ)
    {
        const complexFP * const ptrJ = &(J.toVector()[0]);
        complexFP * const ptrTmpJ = &(tmpJ.toVector()[0]);
        const int n = J.toVector().size();

#pragma omp parallel for
        for (int i = 0; i < n; i++)
            ptrTmpJ[i] = ptrJ[i];
    }

    inline void PSATDTimeStraggered::saveJ()
    {
        assignJ(complexGrid->Jx, tmpJx);
        assignJ(complexGrid->Jy, tmpJy);
        assignJ(complexGrid->Jz, tmpJz);
    }

    inline void PSATDTimeStraggered::updateFields()
    {
        doFourierTransform(RtoC);

        if (pml.get()) getPml()->updateBSplit();
        updateHalfB();

        if (pml.get()) getPml()->updateESplit();
        updateE();

        if (pml.get()) getPml()->updateBSplit();
        updateHalfB();

        saveJ();
        doFourierTransform(CtoR);

        if (pml.get()) getPml()->doSecondStep();

        grid->globalTime += grid->dt;
        complexGrid->globalTime += grid->dt;
    }

    inline void PSATDTimeStraggered::convertFieldsPoissonEquation() {
        doFourierTransform(RtoC);
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;
        double dt = grid->dt / 2;
#pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
            {
                //#pragma omp simd
                for (int k = begin.z; k < end.z; k++)
                {
                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();

                    if (normK == 0) {
                        continue;
                    }

                    K = K / normK;

                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 El = (ComplexFP3)K * dot((ComplexFP3)K, E);

                    complexGrid->Ex(i, j, k) -= El.x;
                    complexGrid->Ey(i, j, k) -= El.y;
                    complexGrid->Ez(i, j, k) -= El.z;
                }
            }
        doFourierTransform(CtoR);
    }

    inline void PSATDTimeStraggered::updateHalfB()
    {
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;
        double dt = grid->dt*0.5;
#pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
            {
                //#pragma omp simd
                for (int k = begin.z; k < end.z; k++)
                {
                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();
                    if (normK == 0) {
                        continue;
                    }
                    K = K / normK;

                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 J(complexGrid->Jx(i, j, k), complexGrid->Jy(i, j, k), complexGrid->Jz(i, j, k)),
                        prevJ(tmpJx(i, j, k), tmpJy(i, j, k), tmpJz(i, j, k));
                    ComplexFP3 crossKE = cross((ComplexFP3)K, E);
                    ComplexFP3 crossKJ = cross((ComplexFP3)K, J - prevJ);

                    FP S = sin(normK*constants::c*dt*0.5), C = cos(normK*constants::c*dt*0.5);
                    complexFP coeff1 = 2 * complexFP::i()*S, coeff2 = complexFP::i() * ((1 - C) / (normK*constants::c));

                    complexGrid->Bx(i, j, k) += -coeff1 * crossKE.x + coeff2 * crossKJ.x;
                    complexGrid->By(i, j, k) += -coeff1 * crossKE.y + coeff2 * crossKJ.y;
                    complexGrid->Bz(i, j, k) += -coeff1 * crossKE.z + coeff2 * crossKJ.z;
                }
            }
    }

    inline void PSATDTimeStraggered::updateE()
    {
        const Int3 begin = updateComplexEAreaBegin;
        const Int3 end = updateComplexEAreaEnd;
        double dt = grid->dt;
#pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
            {
                //#pragma omp simd
                for (int k = begin.z; k < end.z; k++)
                {
                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();
                    if (normK == 0) {
                        complexGrid->Ex(i, j, k) += dt * complexGrid->Jx(i, j, k);
                        complexGrid->Ey(i, j, k) += dt * complexGrid->Jy(i, j, k);
                        complexGrid->Ez(i, j, k) += dt * complexGrid->Jz(i, j, k);
                        continue;
                    }
                    K = K / normK;

                    ComplexFP3 B(complexGrid->Bx(i, j, k), complexGrid->By(i, j, k), complexGrid->Bz(i, j, k));
                    ComplexFP3 J(complexGrid->Jx(i, j, k), complexGrid->Jy(i, j, k), complexGrid->Jz(i, j, k));
                    ComplexFP3 crossKB = cross((ComplexFP3)K, B);
                    ComplexFP3 Jl = (ComplexFP3)K * dot((ComplexFP3)K, J);

                    FP S = sin(normK*constants::c*dt*0.5);
                    complexFP coeff1 = 2 * complexFP::i()*S, coeff2 = 2 * S / (normK*constants::c),
                        coeff3 = coeff2 - dt;

                    complexGrid->Ex(i, j, k) += coeff1 * crossKB.x - coeff2 * J.x + coeff3 * Jl.x;
                    complexGrid->Ey(i, j, k) += coeff1 * crossKB.y - coeff2 * J.y + coeff3 * Jl.y;
                    complexGrid->Ez(i, j, k) += coeff1 * crossKB.z - coeff2 * J.z + coeff3 * Jl.z;
                }
            }
    }


    class PSATD : public SpectralFieldSolver<PSATDGridType>
    {
    public:
        PSATD(PSATDGrid* grid);

        void updateFields();

        virtual void updateEB();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);

        void setTimeStep(FP dt);

        void convertFieldsPoissonEquation();

    private:

        PmlSpectralTimeStraggered<GridTypes::PSATDGridType>* getPml() {
            return (PmlSpectralTimeStraggered<GridTypes::PSATDGridType>*)pml.get();
        }

    };

    inline PSATD::PSATD(PSATDGrid* _grid) :
        SpectralFieldSolver<GridTypes::PSATDGridType>(_grid)
    {
        updateDims();
        updateInternalDims();
    }

    inline void PSATD::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        pml.reset(new PmlPsatd(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    inline void PSATD::setTimeStep(FP dt)
    {
        if (grid->setTimeStep(dt)) {
            complexGrid->setTimeStep(dt);
            if (pml.get()) pml.reset(new PmlPsatd(this, pml->sizePML));
        }
    }

    inline void PSATD::updateFields() {
        // std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        doFourierTransform(RtoC);
        //std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        //std::chrono::milliseconds timeRtoC = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

        //std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
        if (pml.get()) getPml()->updateBSplit();
        updateEB();
        if (pml.get()) getPml()->updateESplit();
        updateEB();
        if (pml.get()) getPml()->updateBSplit();
        //std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
        //std::chrono::milliseconds timeSolver = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);

        //std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
        doFourierTransform(CtoR);
        //std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();
        //std::chrono::milliseconds timeCtoR = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5);

        if (pml.get()) getPml()->doSecondStep();

        grid->globalTime += grid->dt;
        complexGrid->globalTime += grid->dt;

        //std::string strRtoC = "Time RtoC: " + std::to_string(timeRtoC.count()) + "\n";
        //std::string strSolver = "Time PSATD: " + std::to_string(timeSolver.count()) + "\n";
        //std::string strCtoR = "Time CtoR: " + std::to_string(timeCtoR.count()) + "\n";
        //std::cout << strRtoC << strSolver << strCtoR << std::endl;
    }

    inline void PSATD::convertFieldsPoissonEquation() {
        doFourierTransform(RtoC);
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;
        double dt = grid->dt / 2;
#pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
            {
                //#pragma omp simd
                for (int k = begin.z; k < end.z; k++)
                {
                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();

                    if (normK == 0) {
                        continue;
                    }

                    K = K / normK;

                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 El = (ComplexFP3)K * dot((ComplexFP3)K, E);

                    complexGrid->Ex(i, j, k) -= El.x;
                    complexGrid->Ey(i, j, k) -= El.y;
                    complexGrid->Ez(i, j, k) -= El.z;
                }
            }
        doFourierTransform(CtoR);
    }

    inline void PSATD::updateEB()
    {
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;
        double dt = grid->dt / 2;
#pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
            {
                //#pragma omp simd
                for (int k = begin.z; k < end.z; k++)
                {
                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();

                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 B(complexGrid->Bx(i, j, k), complexGrid->By(i, j, k), complexGrid->Bz(i, j, k));
                    ComplexFP3 J(complexGrid->Jx(i, j, k), complexGrid->Jy(i, j, k), complexGrid->Jz(i, j, k));
                    J = complexFP(4 * constants::pi) * J;

                    if (normK == 0) {
                        complexGrid->Ex(i, j, k) += -J.x;
                        complexGrid->Ey(i, j, k) += -J.y;
                        complexGrid->Ez(i, j, k) += -J.z;
                        continue;
                    }

                    K = K / normK;

                    ComplexFP3 kEcross = cross((ComplexFP3)K, E), kBcross = cross((ComplexFP3)K, B),
                        kJcross = cross((ComplexFP3)K, J);
                    ComplexFP3 Jl = (ComplexFP3)K * dot((ComplexFP3)K, J), El = (ComplexFP3)K * dot((ComplexFP3)K, E);

                    FP S = sin(normK*constants::c*dt), C = cos(normK*constants::c*dt);

                    complexFP coef1E = S * complexFP::i(), coef2E = -S / (normK*constants::c),
                        coef3E = S / (normK*constants::c) - dt;

                    complexGrid->Ex(i, j, k) = C * E.x + coef1E * kBcross.x + (1 - C) * El.x + coef2E * J.x + coef3E * Jl.x;
                    complexGrid->Ey(i, j, k) = C * E.y + coef1E * kBcross.y + (1 - C) * El.y + coef2E * J.y + coef3E * Jl.y;
                    complexGrid->Ez(i, j, k) = C * E.z + coef1E * kBcross.z + (1 - C) * El.z + coef2E * J.z + coef3E * Jl.z;

                    complexFP coef1B = -S * complexFP::i(), coef2B = ((1 - C) / (normK*constants::c))*complexFP::i();

                    complexGrid->Bx(i, j, k) = C * B.x + coef1B * kEcross.x + coef2B * kJcross.x;
                    complexGrid->By(i, j, k) = C * B.y + coef1B * kEcross.y + coef2B * kJcross.y;
                    complexGrid->Bz(i, j, k) = C * B.z + coef1B * kEcross.z + coef2B * kJcross.z;
                }
            }
    }
}
