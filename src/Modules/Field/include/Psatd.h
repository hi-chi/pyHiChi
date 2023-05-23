#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "PmlPsatd.h"
#include <omp.h>

namespace pfc {

    template <bool ifPoisson>
    class PSATDT : public SpectralFieldSolver<PSATDGridType>
    {
    public:

        using GridType = PSATDGrid;
        using PmlType = PmlPsatd;
        using PeriodicalFieldGeneratorType = PeriodicalFieldGeneratorPsatd;

        PSATDT(GridType* grid, FP dt);

        void updateFields();

        virtual void updateEB();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void setFieldGenerator(PeriodicalFieldGeneratorType* _generator);

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
    inline PSATDT<ifPoisson>::PSATDT(PSATDT<ifPoisson>::GridType* grid, FP dt) :
        SpectralFieldSolver(grid, dt, 0.0, 0.0, 0.5 * dt)
    {
        updateDims();
        updateInternalDims();
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        pml.reset(new PmlType(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setFieldGenerator(PSATDT<ifPoisson>::PeriodicalFieldGeneratorType* _generator)
    {
        generator.reset(_generator->createInstance(this));
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setTimeStep(FP dt)
    {
        this->dt = dt;
        this->timeShiftJ = 0.5 * dt;
        if (pml) pml.reset(new PmlType(this, pml->sizePml));
        if (generator) generator.reset(generator->createInstance(this));
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::updateFields() {
        doFourierTransform(fourier_transform::Direction::RtoC);

        if (pml) getPml()->updateBSplit();
        updateEB();
        if (pml) getPml()->updateESplit();
        updateEB();
        if (pml) getPml()->updateBSplit();

        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml) getPml()->updateB();
        if (pml) getPml()->updateE();

        globalTime += dt;
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::convertFieldsPoissonEquation() {
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
                            complexGrid->Ex(i, j, k) += -J.x;
                            complexGrid->Ey(i, j, k) += -J.y;
                            complexGrid->Ez(i, j, k) += -J.z;
                            continue;
                        }

                        K = K / normK;

                        ComplexFP3 kEcross = cross((ComplexFP3)K, E), kBcross = cross((ComplexFP3)K, B),
                            kJcross = cross((ComplexFP3)K, J);
                        ComplexFP3 Jl = (ComplexFP3)K * dot((ComplexFP3)K, J), El = (ComplexFP3)K * dot((ComplexFP3)K, E);

                        FP S = sin(normK * constants::c * dt), C = cos(normK * constants::c * dt);

                        complexFP coef1E = S * complexFP::i(), coef2E = -S / (normK * constants::c),
                            coef3E = S / (normK * constants::c) - dt;

                        complexGrid->Ex(i, j, k) = C * E.x + coef1E * kBcross.x + (1 - C) * El.x + coef2E * J.x + coef3E * Jl.x;
                        complexGrid->Ey(i, j, k) = C * E.y + coef1E * kBcross.y + (1 - C) * El.y + coef2E * J.y + coef3E * Jl.y;
                        complexGrid->Ez(i, j, k) = C * E.z + coef1E * kBcross.z + (1 - C) * El.z + coef2E * J.z + coef3E * Jl.z;

                        complexFP coef1B = -S * complexFP::i(), coef2B = ((1 - C) / (normK * constants::c)) * complexFP::i();

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
                            complexGrid->Ex(i, j, k) += -J.x;
                            complexGrid->Ey(i, j, k) += -J.y;
                            complexGrid->Ez(i, j, k) += -J.z;
                            continue;
                        }

                        K = K / normK;

                        ComplexFP3 kEcross = cross((ComplexFP3)K, E), kBcross = cross((ComplexFP3)K, B),
                            kJcross = cross((ComplexFP3)K, J);
                        ComplexFP3 Jl = (ComplexFP3)K * dot((ComplexFP3)K, J), El = (ComplexFP3)K * dot((ComplexFP3)K, E);

                        FP S = sin(normK * constants::c * dt), C = cos(normK * constants::c * dt);

                        complexFP coef1E = S * complexFP::i(), coef2E = -S / (normK * constants::c),
                            coef3E = S / (normK * constants::c) - dt;

                        complexGrid->Ex(i, j, k) = C * (E.x - El.x) + coef1E * kBcross.x + coef2E * (J.x - Jl.x);
                        complexGrid->Ey(i, j, k) = C * (E.y - El.y) + coef1E * kBcross.y + coef2E * (J.y - Jl.y);
                        complexGrid->Ez(i, j, k) = C * (E.z - El.z) + coef1E * kBcross.z + coef2E * (J.z - Jl.z);

                        complexFP coef1B = -S * complexFP::i(), coef2B = ((1 - C) / (normK * constants::c)) * complexFP::i();

                        complexGrid->Bx(i, j, k) = C * B.x + coef1B * kEcross.x + coef2B * kJcross.x;
                        complexGrid->By(i, j, k) = C * B.y + coef1B * kEcross.y + coef2B * kJcross.y;
                        complexGrid->Bz(i, j, k) = C * B.z + coef1B * kEcross.z + coef2B * kJcross.z;
                    }
                }
    }

    typedef PSATDT<true> PSATDPoisson;
    typedef PSATDT<false> PSATD;

}
