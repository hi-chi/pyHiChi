#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "PmlFdtd.h"
#include "FieldBoundaryConditionFdtd.h"
#include "FieldGeneratorFdtd.h"
#include "Vectors.h"

#include <algorithm>

namespace pfc {
    
    class FDTD : public RealFieldSolver<YeeGrid, PmlFdtd, FieldGeneratorFdtd>
    {
    public:

        using GridType = YeeGrid;
        using PmlType = PmlFdtd;
        using FieldGeneratorType = FieldGeneratorFdtd;
        using PeriodicalBoundaryConditionType = PeriodicalBoundaryConditionFdtd;
        using ReflectBoundaryConditionType = ReflectBoundaryConditionFdtd;

        FDTD(GridType* grid, FP dt);

        // constructor for loading
        explicit FDTD(GridType* grid);

        void updateFields();

        void updateHalfB();
        void updateE();

        void setPeriodicalBoundaryConditions();
        void setPeriodicalBoundaryConditions(CoordinateEnum axis);

        void setReflectBoundaryConditions();
        void setReflectBoundaryConditions(CoordinateEnum axis);

        void setTimeStep(FP dt);

        static FP getCourantConditionTimeStep(const FP3& gridSteps) {
            FP tmp = sqrt(1.0 / (gridSteps.x * gridSteps.x) +
                1.0 / (gridSteps.y * gridSteps.y) +
                1.0 / (gridSteps.z * gridSteps.z));
            return 1.0 / (constants::c * tmp);
        }

        FP getCourantConditionTimeStep() const {
            return getCourantConditionTimeStep(this->grid->steps);
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

    private:

        void updateHalfB3D();
        void updateHalfB2D();
        void updateHalfB1D();
        void updateE3D();
        void updateE2D();
        void updateE1D();

        FP3 anisotropyCoeff;
        void setAnisotropy(const FP frequency, int axis);

    };

    inline FDTD::FDTD(GridType* grid, FP dt) :
        RealFieldSolver<GridType, PmlType, FieldGeneratorType>(grid, dt),
        anisotropyCoeff(1, 1, 1)
    {
        if (!isCourantConditionSatisfied(dt)) {
            std::cout
                << "WARNING: FDTD Courant condition is not satisfied. Another time step was setted up"
                << std::endl;
            this->dt = getCourantConditionTimeStep() * 0.5;
        }
    }

    inline FDTD::FDTD(GridType* grid) :
        RealFieldSolver<GridType, PmlType, FieldGeneratorType>(grid)
    {}

    inline void FDTD::setPeriodicalBoundaryConditions()
    {
        for (int d = 0; d < this->grid->dimensionality; d++)
            this->boundaryConditions[d].reset(new PeriodicalBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd, (CoordinateEnum)d));
    }

    inline void FDTD::setPeriodicalBoundaryConditions(CoordinateEnum axis)
    {
        if ((int)axis < this->grid->dimensionality)
            this->boundaryConditions[(int)axis].reset(new PeriodicalBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd, axis));
    }

    inline void FDTD::setReflectBoundaryConditions()
    {
        for (int d = 0; d < this->grid->dimensionality; d++)
            this->boundaryConditions[d].reset(new ReflectBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd));
    }

    inline void FDTD::setReflectBoundaryConditions(CoordinateEnum axis)
    {
        if ((int)axis < this->grid->dimensionality)
            this->boundaryConditions[(int)axis].reset(new ReflectBoundaryConditionType(
                this->grid, this->domainIndexBegin, this->domainIndexEnd, axis));
    }

    inline void FDTD::setTimeStep(FP dt)
    {
        if (isCourantConditionSatisfied(dt)) {
            this->dt = dt;
            if (this->pml) this->pml->dt = dt;
            if (this->generator) this->generator->dt = dt;
        }
        else {
            std::cout
                << "WARNING: FDTD Courant condition is not satisfied. Time step was not changed"
                << std::endl;
        }
    }

    inline void FDTD::setAnisotropy(FP frequency, int axis)
    {
        // We introduce artificial anisotropy, through one axis.
        // For this we upgrade Maxwell equations by coefficients,
        // which computes from major signal frequency. 
        // See more in Juntunen,Tsiboukis - Reduction of Numerical Dispersion in
        // FDTD Method Through Artificial Anisotropy.

        FP3 steps = grid->steps;
        FP WP = constants::pi * 2.0 * constants::c / frequency;
        FP R = WP / steps.norm();
        const FP q = 0.99;  // q - stability coefficient, 0 <= q <= 1
        FP Amax = constants::pi / (3 * R * asin(asin(constants::pi / (R * sqrt(3.0))) / sqrt(3.0)));
        FP Q = Amax - 1;
        FP c1 = 1 - Q / 2;
        int axis0 = axis;
        int axis1 = (axis + 1) % 3;
        int axis2 = (axis + 2) % 3;
        // equivalents of the variables
        // Z1 == Zy, Zz == Zz
        // Zy,Zz - designation from article
        FP Z1 = steps[axis0] / steps[axis1];
        FP Z2 = steps[axis0] / steps[axis2];
        // equivalents of the variables
        // CoeffA == K1, CoeffB == K2, a1 == a, a2 == b
        // K1, K2, a, b - designation from article
        FP CoeffA = constants::pi / (R * sqrt(1 + 1 / (Z1 * Z1) + 1 / (Z2 * Z2)));
        FP a1 = sin(CoeffA / c1) * sin(CoeffA / c1)
            / (Z1 * Z1 * sin(CoeffA / (c1 * Z1)) * sin(CoeffA / (c1 * Z1)));
        FP a2 = sin(CoeffA / c1) * sin(CoeffA / c1)
            / (Z2 * Z2 * sin(CoeffA / (c1 * Z2)) * sin(CoeffA / (c1 * Z2)));
        FP CoeffB = sqrt(1 + a1 * Z1 * Z1 + a2 * Z2 * Z2);
        anisotropyCoeff[axis0] = CoeffB / (CoeffA * q * sqrt(a1 * a2))
            * asin(q * sin(CoeffA / c1) / CoeffB);
        anisotropyCoeff[axis1] = a1 * anisotropyCoeff[axis0];
        anisotropyCoeff[axis2] = a2 * anisotropyCoeff[axis0];
    }

    inline void FDTD::updateFields()
    {
        updateHalfB();
        if (pml) pml->updateB();
        if (generator) generator->generateB(globalTime);  // send current E time
        applyBoundaryConditionsB(globalTime + dt * 0.5);

        updateE();
        if (pml) pml->updateE();
        if (generator) generator->generateE(globalTime + dt * 0.5);  // send current B time
        applyBoundaryConditionsE(globalTime + dt);

        updateHalfB();
        applyBoundaryConditionsB(globalTime + dt);

        globalTime += dt;
    }

    inline void FDTD::updateHalfB()
    {
        if (grid->dimensionality == 3)
            updateHalfB3D();
        else if (grid->dimensionality == 2)
            updateHalfB2D();
        else if (grid->dimensionality == 1)
            updateHalfB1D();
    }

    inline void FDTD::updateHalfB3D()
    {
        const FP cdt = constants::c * dt * (FP)0.5;
        const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
        const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
        const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
        const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);
        const FP coeffZX = cdt / (grid->steps.z * anisotropyCoeff.x);
        const FP coeffZY = cdt / (grid->steps.z * anisotropyCoeff.y);

        const Int3 begin = this->internalIndexBegin;
        const Int3 end = this->internalIndexEnd;

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
            {
                OMP_SIMD()
                for (int k = begin.z; k < end.z; k++)
                {
                    grid->Bx(i, j, k) += coeffZX * (grid->Ey(i, j, k) - grid->Ey(i, j, k - 1)) -
                        coeffYX * (grid->Ez(i, j, k) - grid->Ez(i, j - 1, k));
                    grid->By(i, j, k) += coeffXY * (grid->Ez(i, j, k) - grid->Ez(i - 1, j, k)) -
                        coeffZY * (grid->Ex(i, j, k) - grid->Ex(i, j, k - 1));
                    grid->Bz(i, j, k) += coeffYZ * (grid->Ex(i, j, k) - grid->Ex(i, j - 1, k)) -
                        coeffXZ * (grid->Ey(i, j, k) - grid->Ey(i - 1, j, k));
                }
            }
    }

    inline void FDTD::updateHalfB2D()
    {
        const FP cdt = constants::c * dt * (FP)0.5;
        const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
        const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
        const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
        const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);

        const Int3 begin = this->internalIndexBegin;
        const Int3 end = this->internalIndexEnd;

        OMP_FOR()
        for (int i = begin.x; i < end.x; i++) {
            OMP_SIMD()
            for (int j = begin.y; j < end.y; j++)
            {
                grid->Bx(i, j, 0) += -coeffYX * (grid->Ez(i, j, 0) - grid->Ez(i, j - 1, 0));
                grid->By(i, j, 0) += coeffXY * (grid->Ez(i, j, 0) - grid->Ez(i - 1, j, 0));
                grid->Bz(i, j, 0) += coeffYZ * (grid->Ex(i, j, 0) - grid->Ex(i, j - 1, 0)) -
                    coeffXZ * (grid->Ey(i, j, 0) - grid->Ey(i - 1, j, 0));
            }
        }
    }

    inline void FDTD::updateHalfB1D()
    {
        const FP cdt = constants::c * dt * (FP)0.5;
        const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
        const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);

        const Int3 begin = this->internalIndexBegin;
        const Int3 end = this->internalIndexEnd;

        OMP_FOR()
        for (int i = begin.x; i < end.x; i++) {
            grid->By(i, 0, 0) += coeffXY * (grid->Ez(i, 0, 0) - grid->Ez(i - 1, 0, 0));
            grid->Bz(i, 0, 0) += -coeffXZ * (grid->Ey(i, 0, 0) - grid->Ey(i - 1, 0, 0));
        }
    }

    inline void FDTD::updateE()
    {
        if (grid->dimensionality == 3)
            updateE3D();
        else if (grid->dimensionality == 2)
            updateE2D();
        else if (grid->dimensionality == 1)
            updateE1D();
    }

    inline void FDTD::updateE3D()
    {
        const FP coeffCurrent = -(FP)4 * constants::pi * dt;
        const FP cdt = constants::c * dt;
        const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
        const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
        const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
        const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);
        const FP coeffZX = cdt / (grid->steps.z * anisotropyCoeff.x);
        const FP coeffZY = cdt / (grid->steps.z * anisotropyCoeff.y);

        const Int3 begin = this->internalIndexBegin;
        const Int3 end = this->internalIndexEnd;

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
            {
                OMP_SIMD()
                for (int k = begin.z; k < end.z; k++)
                {
                    grid->Ex(i, j, k) += coeffCurrent * grid->Jx(i, j, k) +
                        coeffYX * (grid->Bz(i, j + 1, k) - grid->Bz(i, j, k)) -
                        coeffZX * (grid->By(i, j, k + 1) - grid->By(i, j, k));
                    grid->Ey(i, j, k) += coeffCurrent * grid->Jy(i, j, k) +
                        coeffZY * (grid->Bx(i, j, k + 1) - grid->Bx(i, j, k)) -
                        coeffXY * (grid->Bz(i + 1, j, k) - grid->Bz(i, j, k));
                    grid->Ez(i, j, k) += coeffCurrent * grid->Jz(i, j, k) +
                        coeffXZ * (grid->By(i + 1, j, k) - grid->By(i, j, k)) -
                        coeffYZ * (grid->Bx(i, j + 1, k) - grid->Bx(i, j, k));
                }
            }
    }

    inline void FDTD::updateE2D()
    {
        const FP coeffCurrent = -(FP)4 * constants::pi * dt;
        const FP cdt = constants::c * dt;
        const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
        const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
        const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
        const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);

        const Int3 begin = this->internalIndexBegin;
        const Int3 end = this->internalIndexEnd;

        OMP_FOR()
        for (int i = begin.x; i < end.x; i++) {
            OMP_SIMD()
            for (int j = begin.y; j < end.y; j++) {
                grid->Ex(i, j, 0) += coeffCurrent * grid->Jx(i, j, 0) +
                    coeffYX * (grid->Bz(i, j + 1, 0) - grid->Bz(i, j, 0));
                grid->Ey(i, j, 0) += coeffCurrent * grid->Jy(i, j, 0) -
                    coeffXY * (grid->Bz(i + 1, j, 0) - grid->Bz(i, j, 0));
                grid->Ez(i, j, 0) += coeffCurrent * grid->Jz(i, j, 0) +
                    coeffXZ * (grid->By(i + 1, j, 0) - grid->By(i, j, 0)) -
                    coeffYZ * (grid->Bx(i, j + 1, 0) - grid->Bx(i, j, 0));
            }
        }
    }

    inline void FDTD::updateE1D()
    {
        const FP coeffCurrent = -(FP)4 * constants::pi * dt;
        const FP cdt = constants::c * dt;
        const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
        const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);

        const Int3 begin = this->internalIndexBegin;
        const Int3 end = this->internalIndexEnd;

        OMP_FOR()
        for (int i = begin.x; i < end.x; i++) {
            grid->Ex(i, 0, 0) += coeffCurrent * grid->Jx(i, 0, 0);
            grid->Ey(i, 0, 0) += coeffCurrent * grid->Jy(i, 0, 0) -
                coeffXY * (grid->Bz(i + 1, 0, 0) - grid->Bz(i, 0, 0));
            grid->Ez(i, 0, 0) += coeffCurrent * grid->Jz(i, 0, 0) +
                coeffXZ * (grid->By(i + 1, 0, 0) - grid->By(i, 0, 0));
        }
    }

    inline void FDTD::save(std::ostream& ostr)
    {
        RealFieldSolver<GridType, PmlType, FieldGeneratorType>::save(ostr);
        ostr.write((char*)&anisotropyCoeff, sizeof(anisotropyCoeff));

        this->saveFieldGenerator(ostr);
        this->savePML(ostr);
        this->saveBoundaryConditions(ostr);
    }

    inline void FDTD::load(std::istream& istr)
    {
        RealFieldSolver<GridType, PmlType, FieldGeneratorType>::load(istr);
        istr.read((char*)&anisotropyCoeff, sizeof(anisotropyCoeff));

        this->loadFieldGenerator(istr);
        this->loadPML(istr);
        this->loadBoundaryConditions(istr);
    }

    inline void FDTD::saveBoundaryConditions(std::ostream& ostr)
    {
        for (int d = 0; d < 3; d++) {
            int isPeriodicalBC = dynamic_cast<PeriodicalBoundaryConditionType*>(this->boundaryConditions[d].get()) ? 1 : 0;
            ostr.write((char*)&isPeriodicalBC, sizeof(isPeriodicalBC));

            int isReflectBC = dynamic_cast<ReflectBoundaryConditionType*>(this->boundaryConditions[d].get()) ? 1 : 0;
            ostr.write((char*)&isReflectBC, sizeof(isReflectBC));

            if (this->boundaryConditions[d])
                this->boundaryConditions[d]->save(ostr);
        }
    }

    inline void FDTD::loadBoundaryConditions(std::istream& istr)
    {
        for (int d = 0; d < 3; d++) {
            int isPeriodicalBC = 0;
            istr.read((char*)&isPeriodicalBC, sizeof(isPeriodicalBC));

            int isReflectBC = 0;
            istr.read((char*)&isReflectBC, sizeof(isReflectBC));

            if (isPeriodicalBC) {
                this->boundaryConditions[d].reset(new PeriodicalBoundaryConditionType(
                    this->grid, this->domainIndexBegin, this->domainIndexEnd));
                this->boundaryConditions[d]->load(istr);
            }
            else if (isReflectBC) {
                this->boundaryConditions[d].reset(new ReflectBoundaryConditionType(
                    this->grid, this->domainIndexBegin, this->domainIndexEnd));
                this->boundaryConditions[d]->load(istr);
            }
        }
    }
}
