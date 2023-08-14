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
    class PSATDT : public SpectralFieldSolver<PSATDGrid, PmlPsatdTimeStaggered<PSATDGrid>,
        FieldGeneratorSpectral<PSATDGrid>>
    {
    public:

        using GridType = PSATDGrid;
        using PmlType = PmlPsatdTimeStaggered<PSATDGrid>;
        using FieldGeneratorType = FieldGeneratorSpectral<PSATDGrid>;
        using PeriodicalBoundaryConditionType = PeriodicalBoundaryConditionSpectral<PSATDGrid>;

        PSATDT(GridType* grid, FP dt);

        // constructor for loading
        explicit PSATDT(GridType* grid);

        void updateFields();

        void updateEB();

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
            return false;
        }

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        void saveBoundaryConditions(std::ostream& ostr);
        void loadBoundaryConditions(std::istream& istr);

    };

    typedef PSATDT<true> PSATDPoisson;
    typedef PSATDT<false> PSATD;

    template <bool ifPoisson>
    inline PSATDT<ifPoisson>::PSATDT(GridType* grid, FP dt) :
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>(grid, dt)
    {}

    template <bool ifPoisson>
    inline PSATDT<ifPoisson>::PSATDT(GridType* grid) :
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>(grid)
    {}

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setPeriodicalBoundaryConditions()
    {
        for (int d = 0; d < this->grid->dimensionality; d++)
            this->boundaryConditions[d].reset(new PeriodicalBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd, (CoordinateEnum)d));
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setPeriodicalBoundaryConditions(CoordinateEnum axis)
    {
        if ((int)axis < this->grid->dimensionality)
            this->boundaryConditions[(int)axis].reset(new PeriodicalBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd, axis));
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setTimeStep(FP dt)
    {
        this->dt = dt;
        if (this->pml) this->pml->dt = dt;
        if (this->generator) this->generator->dt = dt;
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::updateFields()
    {
        // TODO: consider boundary conditions and generator

        doFourierTransform(fourier_transform::Direction::RtoC);

        if (pml) pml->updateBSplit();
        updateEB();
        // if (generator) generator->generateB(globalTime + dt * 0.5);
        // if (generator) generator->generateE(globalTime + dt * 0.5);
        // applyBoundaryConditionsB(globalTime + dt * 0.5);
        // applyBoundaryConditionsE(globalTime + dt * 0.5);

        if (pml) pml->updateESplit();
        updateEB();
        // if (generator) generator->generateB(globalTime + dt);
        // if (generator) generator->generateE(globalTime + dt);
        // applyBoundaryConditionsB(globalTime + dt);
        // applyBoundaryConditionsE(globalTime + dt);

        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml) pml->updateB();
        if (pml) pml->updateE();

        globalTime += dt;
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::convertFieldsPoissonEquation()
    {
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
    inline void PSATDT<ifPoisson>::updateEB()
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

    // provides k \cdot E = 0 always (k \cdot J = 0 too)
    template <>
    inline void PSATDT<true>::updateEB()
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

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::save(std::ostream& ostr)
    {
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>::save(ostr);

        this->saveFieldGenerator(ostr);
        this->savePML(ostr);
        this->saveBoundaryConditions(ostr);
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::load(std::istream& istr)
    {
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>::load(istr);

        this->loadFieldGenerator(istr);
        this->loadPML(istr);
        this->loadBoundaryConditions(istr);
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::saveBoundaryConditions(std::ostream& ostr)
    {
        for (int d = 0; d < 3; d++) {
            int isPeriodicalBC = dynamic_cast<PeriodicalBoundaryConditionType*>(this->boundaryConditions[d].get()) ? 1 : 0;
            ostr.write((char*)&isPeriodicalBC, sizeof(isPeriodicalBC));

            if (this->boundaryConditions[d])
                this->boundaryConditions[d]->save(ostr);
        }
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::loadBoundaryConditions(std::istream& istr)
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
