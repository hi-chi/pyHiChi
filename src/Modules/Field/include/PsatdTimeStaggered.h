#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "PmlPsatd.h"
#include "FieldBoundaryConditionSpectral.h"
#include "FieldGeneratorSpectral.h"

namespace pfc {

    template <bool ifPoisson>
    class PSATDTimeStaggeredT : public SpectralFieldSolver<PSATDTimeStaggeredGrid,
        PmlPsatdTimeStaggered<PSATDTimeStaggeredGrid>, FieldGeneratorSpectral<PSATDTimeStaggeredGrid>>
    {
    public:

        using GridType = PSATDTimeStaggeredGrid;
        using PmlType = PmlPsatdTimeStaggered<PSATDTimeStaggeredGrid>;
        using FieldGeneratorType = FieldGeneratorSpectral<PSATDTimeStaggeredGrid>;
        using PeriodicalBoundaryConditionType = PeriodicalBoundaryConditionSpectral<PSATDTimeStaggeredGrid>;

        PSATDTimeStaggeredT(GridType* grid);  // use when load
        PSATDTimeStaggeredT(GridType* grid, FP dt);

        void updateFields();

        void updateHalfB();
        void updateE();

        void convertFieldsPoissonEquation();

        void setTimeStep(FP dt);

        static FP getCourantConditionTimeStep(const FP3& gridSteps) {
            // PSATD is free of the Courant condition
            // but PIC-code requires the following condition
            FP tmp = sqrt(1.0 / (gridSteps.x * gridSteps.x) +
                1.0 / (gridSteps.y * gridSteps.y) +
                1.0 / (gridSteps.z * gridSteps.z));
            return 1.0 / (constants::c * tmp);
        }
        FP getCourantConditionTimeStep() const {
            return getCourantConditionTimeStep(grid->steps);
        }

        bool isCourantConditionSatisfied(FP dt) const {
            return true;
        }

        static bool isTimeStaggered() {
            return true;
        }

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        ScalarField<complexFP> tmpJx, tmpJy, tmpJz;

    protected:

        void saveJ();
        void assignJ(SpectralScalarField<FP, complexFP>& J, ScalarField<complexFP>& tmpJ);
    };

    template <bool ifPoisson>
    inline PSATDTimeStaggeredT<ifPoisson>::PSATDTimeStaggeredT(GridType* grid) :
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>(grid),
        tmpJx(this->complexGrid->sizeStorage),
        tmpJy(this->complexGrid->sizeStorage),
        tmpJz(this->complexGrid->sizeStorage)
    {}

    template <bool ifPoisson>
    inline PSATDTimeStaggeredT<ifPoisson>::PSATDTimeStaggeredT(GridType* grid, FP dt) :
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>(grid, dt),
        tmpJx(this->complexGrid->sizeStorage),
        tmpJy(this->complexGrid->sizeStorage),
        tmpJz(this->complexGrid->sizeStorage)
    {}

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::setTimeStep(FP dt)
    {
        this->dt = dt;
        resetPML();
        resetFieldGenerator();
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::assignJ(
        SpectralScalarField<FP, complexFP>& J, ScalarField<complexFP>& tmpJ)
    {
        const complexFP* const ptrJ = J.getData();
        complexFP* const ptrTmpJ = tmpJ.getData();
        const int n = J.getSize().volume();

        OMP_FOR()
        for (int i = 0; i < n; i++)
            ptrTmpJ[i] = ptrJ[i];
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::saveJ()
    {
        assignJ(complexGrid->Jx, tmpJx);
        assignJ(complexGrid->Jy, tmpJy);
        assignJ(complexGrid->Jz, tmpJz);
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::updateFields()
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

        if (pml) pml->updateBSplit();
        updateHalfB();
        // applyBoundaryConditionsB(globalTime + dt);

        saveJ();
        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml) pml->updateB();
        if (pml) pml->updateE();

        globalTime += dt;
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::convertFieldsPoissonEquation() {
        doFourierTransform(fourier_transform::Direction::RtoC);
        
        const Int3 begin = this->complexDomainIndexBegin;
        const Int3 end = this->complexDomainIndexEnd;

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
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
        doFourierTransform(fourier_transform::Direction::CtoR);
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::updateHalfB()
    {
        const Int3 begin = this->complexDomainIndexBegin;
        const Int3 end = this->complexDomainIndexEnd;

        FP dt = 0.5 * this->dt;

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
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

                    FP S = sin(normK * constants::c * dt * 0.5), C = cos(normK * constants::c * dt * 0.5);
                    complexFP coeff1 = 2 * complexFP::i() * S, coeff2 = complexFP::i() * ((1 - C) / (normK * constants::c));

                    complexGrid->Bx(i, j, k) += -coeff1 * crossKE.x + coeff2 * crossKJ.x;
                    complexGrid->By(i, j, k) += -coeff1 * crossKE.y + coeff2 * crossKJ.y;
                    complexGrid->Bz(i, j, k) += -coeff1 * crossKE.z + coeff2 * crossKJ.z;
                }
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::updateE()
    {
        const Int3 begin = this->complexDomainIndexBegin;
        const Int3 end = this->complexDomainIndexEnd;

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                { 
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

                    ComplexFP3 B(complexGrid->Bx(i, j, k), complexGrid->By(i, j, k), complexGrid->Bz(i, j, k));
                    ComplexFP3 crossKB = cross((ComplexFP3)K, B);
                    ComplexFP3 Jl = (ComplexFP3)K * dot((ComplexFP3)K, J);

                    FP S = sin(normK * constants::c * dt * 0.5);
                    complexFP coeff1 = 2 * complexFP::i() * S, coeff2 = 2 * S / (normK * constants::c),
                        coeff3 = coeff2 - dt;

                    complexGrid->Ex(i, j, k) += coeff1 * crossKB.x - coeff2 * J.x + coeff3 * Jl.x;
                    complexGrid->Ey(i, j, k) += coeff1 * crossKB.y - coeff2 * J.y + coeff3 * Jl.y;
                    complexGrid->Ez(i, j, k) += coeff1 * crossKB.z - coeff2 * J.z + coeff3 * Jl.z;
                }
    }

    // provides k \cdot E = 0 always (k \cdot J = 0 too)
    template <>
    inline void PSATDTimeStaggeredT<true>::updateE()
    {
        const Int3 begin = this->complexDomainIndexBegin;
        const Int3 end = this->complexDomainIndexEnd;

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
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

                    ComplexFP3 E(complexGrid->Ex(i, j, k), complexGrid->Ey(i, j, k), complexGrid->Ez(i, j, k));
                    ComplexFP3 B(complexGrid->Bx(i, j, k), complexGrid->By(i, j, k), complexGrid->Bz(i, j, k));
                    ComplexFP3 crossKB = cross((ComplexFP3)K, B);
                    ComplexFP3 El = (ComplexFP3)K * dot((ComplexFP3)K, E);
                    ComplexFP3 Jl = (ComplexFP3)K * dot((ComplexFP3)K, J);

                    FP S = sin(normK * constants::c * dt * 0.5);
                    complexFP coeff1 = 2 * complexFP::i() * S, coeff2 = 2 * S / (normK * constants::c),
                        coeff3 = coeff2 - dt;
                    complexGrid->Ex(i, j, k) += -El.x + coeff1 * crossKB.x - coeff2 * (J.x - Jl.x);
                    complexGrid->Ey(i, j, k) += -El.y + coeff1 * crossKB.y - coeff2 * (J.y - Jl.y);
                    complexGrid->Ez(i, j, k) += -El.z + coeff1 * crossKB.z - coeff2 * (J.z - Jl.z);
                }
    }

    template<bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::save(std::ostream& ostr)
    {
        FieldSolver<GridType, PmlType, FieldGeneratorType>::save(ostr);
        tmpJx.save(ostr);
        tmpJy.save(ostr);
        tmpJz.save(ostr);
    }

    template<bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::load(std::istream& istr)
    {
        FieldSolver<GridType, PmlType, FieldGeneratorType>::load(istr);
        tmpJx.load(istr);
        tmpJy.load(istr);
        tmpJz.load(istr);
    }

    typedef PSATDTimeStaggeredT<true> PSATDTimeStaggeredPoisson;
    typedef PSATDTimeStaggeredT<false> PSATDTimeStaggered;

}
