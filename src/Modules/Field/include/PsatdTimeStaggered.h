#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "PmlPsatd.h"
#include "FieldBoundaryConditionSpectral.h"
#include "FieldGeneratorSpectral.h"

namespace pfc {

    namespace psatd_time_staggered {
        struct SchemeParams {
            using GridType = PSATDTimeStaggeredGrid;
            using PmlType = PmlPsatdTimeStaggered<PSATDTimeStaggeredGrid>;
            using FieldGeneratorType = FieldGeneratorSpectral<PSATDTimeStaggeredGrid>;
            using PeriodicalBoundaryConditionType = PeriodicalBoundaryConditionSpectral<PSATDTimeStaggeredGrid>;
        };
    }

    template <bool ifPoisson>
    class PSATDTimeStaggeredT : public SpectralFieldSolver<psatd_time_staggered::SchemeParams>
    {
    public:

        using GridType = psatd_time_staggered::SchemeParams::GridType;
        using PmlType = psatd_time_staggered::SchemeParams::PmlType;
        using FieldGeneratorType = psatd_time_staggered::SchemeParams::FieldGeneratorType;
        using PeriodicalBoundaryConditionType = psatd_time_staggered::SchemeParams::PeriodicalBoundaryConditionType;

        PSATDTimeStaggeredT(GridType* grid, FP dt);

        // constructor for loading
        explicit PSATDTimeStaggeredT(GridType* grid);

        void updateFields();

        void updateHalfB();
        void updateE();

        void convertFieldsPoissonEquation();

        void setPeriodicalBoundaryConditions();
        void setPeriodicalBoundaryConditions(CoordinateEnum axis);

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

        void saveBoundaryConditions(std::ostream& ostr);
        void loadBoundaryConditions(std::istream& istr);

        ScalarField<complexFP> tmpJx, tmpJy, tmpJz;

    protected:

        void saveJ();
        void assignJ(SpectralScalarField<FP, complexFP>& J, ScalarField<complexFP>& tmpJ);
    };

    typedef PSATDTimeStaggeredT<true> PSATDTimeStaggeredPoisson;
    typedef PSATDTimeStaggeredT<false> PSATDTimeStaggered;

    template <bool ifPoisson>
    inline PSATDTimeStaggeredT<ifPoisson>::PSATDTimeStaggeredT(GridType* grid, FP dt) :
        SpectralFieldSolver<psatd_time_staggered::SchemeParams>(grid, dt),
        tmpJx(this->complexGrid->numCells),
        tmpJy(this->complexGrid->numCells),
        tmpJz(this->complexGrid->numCells)
    {}

    template <bool ifPoisson>
    inline PSATDTimeStaggeredT<ifPoisson>::PSATDTimeStaggeredT(GridType* grid) :
        SpectralFieldSolver<psatd_time_staggered::SchemeParams>(grid),
        tmpJx(this->complexGrid->numCells),
        tmpJy(this->complexGrid->numCells),
        tmpJz(this->complexGrid->numCells)
    {}

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::setPeriodicalBoundaryConditions()
    {
        for (int d = 0; d < this->grid->dimensionality; d++)
            this->boundaryConditions[d].reset(new PeriodicalBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd, (CoordinateEnum)d));
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::setPeriodicalBoundaryConditions(CoordinateEnum axis)
    {
        if ((int)axis < this->grid->dimensionality)
            this->boundaryConditions[(int)axis].reset(new PeriodicalBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd, axis));
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::setTimeStep(FP dt)
    {
        this->dt = dt;
        if (this->pml) this->pml->dt = dt;
        if (this->generator) this->generator->dt = dt;
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::assignJ(
        SpectralScalarField<FP, complexFP>& J, ScalarField<complexFP>& tmpJ)
    {
        const Int3 begin = Int3(0, 0, 0);
        const Int3 end = J.getSize();

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                    tmpJ(i, j, k) = J(i, j, k);
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

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::save(std::ostream& ostr)
    {
        SpectralFieldSolver<psatd_time_staggered::SchemeParams>::save(ostr);

        tmpJx.save(ostr);
        tmpJy.save(ostr);
        tmpJz.save(ostr);

        this->saveFieldGenerator(ostr);
        this->savePML(ostr);
        this->saveBoundaryConditions(ostr);
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::load(std::istream& istr)
    {
        SpectralFieldSolver<psatd_time_staggered::SchemeParams>::load(istr);
        
        tmpJx.load(istr);
        tmpJy.load(istr);
        tmpJz.load(istr);

        this->loadFieldGenerator(istr);
        this->loadPML(istr);
        this->loadBoundaryConditions(istr);
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::saveBoundaryConditions(std::ostream& ostr)
    {
        for (int d = 0; d < 3; d++) {
            int isPeriodicalBC = dynamic_cast<PeriodicalBoundaryConditionType*>(this->boundaryConditions[d].get()) ? 1 : 0;
            ostr.write((char*)&isPeriodicalBC, sizeof(isPeriodicalBC));

            if (this->boundaryConditions[d])
                this->boundaryConditions[d]->save(ostr);
        }
    }

    template <bool ifPoisson>
    inline void PSATDTimeStaggeredT<ifPoisson>::loadBoundaryConditions(std::istream& istr)
    {
        for (int d = 0; d < 3; d++) {
            int isPeriodicalBC = 0;
            istr.read((char*)&isPeriodicalBC, sizeof(isPeriodicalBC));

            if (isPeriodicalBC) {
                this->boundaryConditions[d].reset(new PeriodicalBoundaryConditionType(
                    this->grid, this->domainIndexBegin, this->domainIndexEnd));
                this->boundaryConditions[d]->load(istr);
            }
        }
    }

}
