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

        PSTD(GridType* grid, FP dt);

        // constructor for loading
        explicit PSTD(GridType* grid);

        void updateFields();

        void updateHalfB();
        void updateE();

        void setPeriodicalBoundaryConditions();
        void setPeriodicalBoundaryConditions(CoordinateEnum axis);

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

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        void saveBoundaryConditions(std::ostream& ostr);
        void loadBoundaryConditions(std::istream& istr);
    };

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

    inline PSTD::PSTD(GridType* grid) :
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>(grid)
    {}

    inline void PSTD::setPeriodicalBoundaryConditions()
    {
        for (int d = 0; d < this->grid->dimensionality; d++)
            this->boundaryConditions[d].reset(new PeriodicalBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd, (CoordinateEnum)d));
    }

    inline void PSTD::setPeriodicalBoundaryConditions(CoordinateEnum axis)
    {
        if ((int)axis < this->grid->dimensionality)
            this->boundaryConditions[(int)axis].reset(new PeriodicalBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd, axis));
    }

    inline void PSTD::setTimeStep(FP dt)
    {
        if (isCourantConditionSatisfied(dt)) {
            this->dt = dt;
            if (this->pml) this->pml->dt = dt;
            if (this->generator) this->generator->dt = dt;
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

    inline void PSTD::save(std::ostream& ostr)
    {
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>::save(ostr);

        this->saveFieldGenerator(ostr);
        this->savePML(ostr);
        this->saveBoundaryConditions(ostr);
    }

    inline void PSTD::load(std::istream& istr)
    {
        SpectralFieldSolver<GridType, PmlType, FieldGeneratorType>::load(istr);

        this->loadFieldGenerator(istr);
        this->loadPML(istr);
        this->loadBoundaryConditions(istr);
    }

    inline void PSTD::saveBoundaryConditions(std::ostream& ostr)
    {
        for (int d = 0; d < 3; d++) {
            int isPeriodicalBC = dynamic_cast<PeriodicalBoundaryConditionType*>(this->boundaryConditions[d].get()) ? 1 : 0;
            ostr.write((char*)&isPeriodicalBC, sizeof(isPeriodicalBC));

            if (this->boundaryConditions[d])
                this->boundaryConditions[d]->save(ostr);
        }
    }

    inline void PSTD::loadBoundaryConditions(std::istream& istr)
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
