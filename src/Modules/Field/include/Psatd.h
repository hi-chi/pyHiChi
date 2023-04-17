#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "PmlPsatd.h"
//#include <chrono>
#include <omp.h>

namespace pfc {

    template <bool ifPoisson>
    class PSATDTimeStraggeredT : public SpectralFieldSolver<PSATDTimeStraggeredGridType>
    {

    public:
        PSATDTimeStraggeredT(PSATDTimeStraggeredGrid * grid, FP dt);

        void updateFields();

        void updateHalfB();
        void updateE();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);

        void setTimeStep(FP dt);

        void convertFieldsPoissonEquation();

        ScalarField<complexFP> tmpJx, tmpJy, tmpJz;

        bool ifCourantConditionSatisfied(FP dt) {
            return true;
        }

    protected:

        PmlSpectral<GridTypes::PSATDTimeStraggeredGridType>* getPml() {
            return (PmlSpectral<GridTypes::PSATDTimeStraggeredGridType>*)pml.get();
        }

        void saveJ();
        void assignJ(ScalarField<complexFP>& J, ScalarField<complexFP>& tmpJ);
    };

    template <bool ifPoisson>
    inline PSATDTimeStraggeredT<ifPoisson>::PSATDTimeStraggeredT(PSATDTimeStraggeredGrid* _grid, FP dt) :
        SpectralFieldSolver<GridTypes::PSATDTimeStraggeredGridType>(_grid, dt, 0.0, 0.5*dt, 0.5*dt),
        tmpJx(complexGrid->sizeStorage),
        tmpJy(complexGrid->sizeStorage),
        tmpJz(complexGrid->sizeStorage)
    {
        updateDims();
        updateInternalDims();
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        pml.reset(new PmlPsatdTimeStraggered(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::setTimeStep(FP dt)
    {
        this->dt = dt;
        this->timeShiftB = 0.5*dt;
        this->timeShiftJ = 0.5*dt;
        if (pml.get()) pml.reset(new PmlPsatdTimeStraggered(this, pml->sizePML));
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::assignJ(ScalarField<complexFP>& J, ScalarField<complexFP>& tmpJ)
    {
        const complexFP * const ptrJ = J.getData();
        complexFP * const ptrTmpJ = tmpJ.getData();
        const int n = J.getSize().volume();

        OMP_FOR()
        for (int i = 0; i < n; i++)
            ptrTmpJ[i] = ptrJ[i];
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::saveJ()
    {
        assignJ(complexGrid->Jx, tmpJx);
        assignJ(complexGrid->Jy, tmpJy);
        assignJ(complexGrid->Jz, tmpJz);
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::updateFields()
    {
        doFourierTransform(fourier_transform::Direction::RtoC);

        if (pml.get()) getPml()->updateBSplit();
        updateHalfB();

        if (pml.get()) getPml()->updateESplit();
        updateE();

        if (pml.get()) getPml()->updateBSplit();
        updateHalfB();

        saveJ();
        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml.get()) getPml()->doSecondStep();

        globalTime += dt;
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::convertFieldsPoissonEquation() {
        doFourierTransform(fourier_transform::Direction::RtoC);
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;
        double dt = this->dt * 0.5;
        OMP_FOR_COLLAPSE()
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
        doFourierTransform(fourier_transform::Direction::CtoR);
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::updateHalfB()
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
                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();
                    if (normK == 0) {
                        continue;
                    }
                    K = K / normK;

                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 J(complexGrid->Jx(i, j, k), complexGrid->Jy(i, j, k), complexGrid->Jz(i, j, k)),
                        prevJ(tmpJx(i, j, k), tmpJy(i, j, k), tmpJz(i, j, k));
                    J = complexFP(4 * constants::pi) * J;
                    prevJ = complexFP(4 * constants::pi) * prevJ;
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

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::updateE()
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
                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 B(complexGrid->Bx(i, j, k), complexGrid->By(i, j, k), complexGrid->Bz(i, j, k));
                    ComplexFP3 J(complexGrid->Jx(i, j, k), complexGrid->Jy(i, j, k), complexGrid->Jz(i, j, k));
                    J = complexFP(4 * constants::pi) * J;

                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();
                    if (normK == 0) {
                        complexGrid->Ex(i, j, k) += -J.x * dt;
                        complexGrid->Ey(i, j, k) += -J.y * dt;
                        complexGrid->Ez(i, j, k) += -J.z * dt;
                        continue;
                    }
                    K = K / normK;

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

    // provides k \cdot E = 0 always (k \cdot J = 0 too)
    template <>
    inline void PSATDTimeStraggeredT<true>::updateE()
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
                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 B(complexGrid->Bx(i, j, k), complexGrid->By(i, j, k), complexGrid->Bz(i, j, k));
                    ComplexFP3 J(complexGrid->Jx(i, j, k), complexGrid->Jy(i, j, k), complexGrid->Jz(i, j, k));
                    J = complexFP(4 * constants::pi) * J;

                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();
                    if (normK == 0) {
                        complexGrid->Ex(i, j, k) += -J.x * dt;
                        complexGrid->Ey(i, j, k) += -J.y * dt;
                        complexGrid->Ez(i, j, k) += -J.z * dt;
                        continue;
                    }
                    K = K / normK;

                    ComplexFP3 crossKB = cross((ComplexFP3)K, B);
                    ComplexFP3 El = (ComplexFP3)K * dot((ComplexFP3)K, E);
                    ComplexFP3 Jl = (ComplexFP3)K * dot((ComplexFP3)K, J);

                    FP S = sin(normK*constants::c*dt*0.5);
                    complexFP coeff1 = 2 * complexFP::i()*S, coeff2 = 2 * S / (normK*constants::c),
                        coeff3 = coeff2 - dt;
                    complexGrid->Ex(i, j, k) += -El.x + coeff1 * crossKB.x - coeff2 * (J.x - Jl.x);
                    complexGrid->Ey(i, j, k) += -El.y + coeff1 * crossKB.y - coeff2 * (J.y - Jl.y);
                    complexGrid->Ez(i, j, k) += -El.z + coeff1 * crossKB.z - coeff2 * (J.z - Jl.z);
                }
            }
    }

    template <bool ifPoisson>
    class PSATDT : public SpectralFieldSolver<PSATDGridType>
    {
    public:
        PSATDT(PSATDGrid* grid, FP dt);

        void updateFields();

        virtual void updateEB();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);

        void setTimeStep(FP dt);

        void convertFieldsPoissonEquation();

        bool ifCourantConditionSatisfied(FP dt) {
            return true;
        }

    private:

        PmlSpectralTimeStraggered<GridTypes::PSATDGridType>* getPml() {
            return (PmlSpectralTimeStraggered<GridTypes::PSATDGridType>*)pml.get();
        }

    };

    template <bool ifPoisson>
    inline PSATDT<ifPoisson>::PSATDT(PSATDGrid* _grid, FP dt) :
        SpectralFieldSolver<GridTypes::PSATDGridType>(_grid, dt, 0.0, 0.0, 0.5*dt)
    {
        updateDims();
        updateInternalDims();
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        pml.reset(new PmlPsatd(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setTimeStep(FP dt)
    {
        this->dt = dt;
        this->timeShiftJ = 0.5*dt;
        if (pml.get()) pml.reset(new PmlPsatd(this, pml->sizePML));
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::updateFields() {
        doFourierTransform(fourier_transform::Direction::RtoC);
        
        if (pml.get()) getPml()->updateBSplit();
        updateEB();
        if (pml.get()) getPml()->updateESplit();
        updateEB();
        if (pml.get()) getPml()->updateBSplit();

        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml.get()) getPml()->doSecondStep();

        globalTime += dt;
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::convertFieldsPoissonEquation() {
        doFourierTransform(fourier_transform::Direction::RtoC);
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;
        double dt = this->dt *0.5;
        OMP_FOR_COLLAPSE()
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
        doFourierTransform(fourier_transform::Direction::CtoR);
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::updateEB()
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
                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();

                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 B(complexGrid->Bx(i, j, k), complexGrid->By(i, j, k), complexGrid->Bz(i, j, k));
                    ComplexFP3 J(complexGrid->Jx(i, j, k), complexGrid->Jy(i, j, k), complexGrid->Jz(i, j, k));
                    J = complexFP(4 * constants::pi) * J;

                    if (normK == 0) {
                        complexGrid->Ex(i, j, k) += -J.x * dt;
                        complexGrid->Ey(i, j, k) += -J.y * dt;
                        complexGrid->Ez(i, j, k) += -J.z * dt;
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

    // provides k \cdot E = 0 always (k \cdot J = 0 too)
    template <>
    inline void PSATDT<true>::updateEB()
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
                    FP3 K = getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();

                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 B(complexGrid->Bx(i, j, k), complexGrid->By(i, j, k), complexGrid->Bz(i, j, k));
                    ComplexFP3 J(complexGrid->Jx(i, j, k), complexGrid->Jy(i, j, k), complexGrid->Jz(i, j, k));
                    J = complexFP(4 * constants::pi) * J;

                    if (normK == 0) {
                        complexGrid->Ex(i, j, k) += -J.x * dt;
                        complexGrid->Ey(i, j, k) += -J.y * dt;
                        complexGrid->Ez(i, j, k) += -J.z * dt;
                        continue;
                    }

                    K = K / normK;

                    ComplexFP3 kEcross = cross((ComplexFP3)K, E), kBcross = cross((ComplexFP3)K, B),
                        kJcross = cross((ComplexFP3)K, J);
                    ComplexFP3 Jl = (ComplexFP3)K * dot((ComplexFP3)K, J), El = (ComplexFP3)K * dot((ComplexFP3)K, E);

                    FP S = sin(normK*constants::c*dt), C = cos(normK*constants::c*dt);

                    complexFP coef1E = S * complexFP::i(), coef2E = -S / (normK*constants::c),
                        coef3E = S / (normK*constants::c) - dt;

                    complexGrid->Ex(i, j, k) = C * (E.x - El.x) + coef1E * kBcross.x + coef2E * (J.x - Jl.x);
                    complexGrid->Ey(i, j, k) = C * (E.y - El.y) + coef1E * kBcross.y + coef2E * (J.y - Jl.y);
                    complexGrid->Ez(i, j, k) = C * (E.z - El.z) + coef1E * kBcross.z + coef2E * (J.z - Jl.z);

                    complexFP coef1B = -S * complexFP::i(), coef2B = ((1 - C) / (normK*constants::c))*complexFP::i();

                    complexGrid->Bx(i, j, k) = C * B.x + coef1B * kEcross.x + coef2B * kJcross.x;
                    complexGrid->By(i, j, k) = C * B.y + coef1B * kEcross.y + coef2B * kJcross.y;
                    complexGrid->Bz(i, j, k) = C * B.z + coef1B * kEcross.z + coef2B * kJcross.z;
                }
            }
    }

    typedef PSATDT<true> PSATDPoisson;
    typedef PSATDT<false> PSATD;
    typedef PSATDTimeStraggeredT<true> PSATDTimeStraggeredPoisson;
    typedef PSATDTimeStraggeredT<false> PSATDTimeStraggered;

}
