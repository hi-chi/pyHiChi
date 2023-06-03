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
    
    class FDTD : public RealFieldSolver<YeeGridType>
    {
    public:

        using GridType = YeeGrid;
        using PmlType = PmlFdtd;
        using FieldGeneratorType = FieldGeneratorFdtd;
        using PeriodicalBoundaryConditionType = PeriodicalBoundaryConditionFdtd;
        using ReflectBoundaryConditionType = ReflectBoundaryConditionFdtd;

        FDTD(GridType* grid, FP dt);

        void updateFields();

        void updateHalfB();
        void updateE();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void setBoundaryCondition(
            FieldBoundaryCondition<GridTypes::YeeGridType>* _boundaryCondition);
        
        void setFieldGenerator(
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FieldGeneratorType::FunctionType bxFunc, FieldGeneratorType::FunctionType byFunc,
            FieldGeneratorType::FunctionType bzFunc, FieldGeneratorType::FunctionType exFunc,
            FieldGeneratorType::FunctionType eyFunc, FieldGeneratorType::FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));
        void setFieldGenerator(
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            /* first index is index of edge (x, y, z),
            second index is index of field component (ex, ey, ez or bx, by, bz) */
            const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftBFunc,
            const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightBFunc,
            const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftEFunc,
            const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightEFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));

        void setTimeStep(FP dt);

        FP getCourantCondition() const {
            FP tmp = sqrt(1.0 / (grid->steps.x*grid->steps.x) +
                1.0 / (grid->steps.y*grid->steps.y) +
                1.0 / (grid->steps.z*grid->steps.z));
            return 1.0 / (constants::c * tmp);
        }

        bool ifCourantConditionSatisfied(FP dt) const {
            return dt < getCourantCondition();
        }

        void updateDims() {
            this->updateEAreaBegin = Int3(0, 0, 0);
            this->updateEAreaEnd = grid->numCells -
                grid->correctNumCellsAccordingToDim(Int3(1, 1, 1));
            this->updateBAreaBegin =
                grid->correctNumCellsAccordingToDim(Int3(1, 1, 1));
            this->updateBAreaEnd = grid->numCells;
        }

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

    inline FDTD::FDTD(FDTD::GridType* grid, FP dt) :
        RealFieldSolver(grid, dt)
    {
        if (!ifCourantConditionSatisfied(dt)) {
            std::cout
                << "WARNING: FDTD Courant condition is not satisfied. Another time step was setted up"
                << std::endl;
            this->dt = getCourantCondition() * 0.5;
        }
        updateDims();
        updateInternalDims();
        anisotropyCoeff = FP3(1, 1, 1);
    }

    inline void FDTD::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        pml.reset(new PmlFdtd(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    inline void FDTD::setBoundaryCondition(
        FieldBoundaryCondition<GridTypes::YeeGridType>* _boundaryCondition)
    {
        boundaryCondition.reset(_boundaryCondition->createInstance(this));
    }

    inline void FDTD::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        FieldGeneratorType::FunctionType bxFunc, FieldGeneratorType::FunctionType byFunc,
        FieldGeneratorType::FunctionType bzFunc, FieldGeneratorType::FunctionType exFunc,
        FieldGeneratorType::FunctionType eyFunc, FieldGeneratorType::FunctionType ezFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        //generator.reset(new FDTD::FieldGeneratorType(
        //    this, leftGenIndex, rightGenIndex,
        //    bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
        //    isLeftBorderEnabled, isRightBorderEnabled)
        //);
        generator.reset(new FDTD::FieldGeneratorType(this));
    }

    inline void FDTD::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftBFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightBFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftEFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightEFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        //generator.reset(new FDTD::FieldGeneratorType(
        //    this, leftGenIndex, rightGenIndex,
        //    leftBFunc, rightBFunc, leftEFunc, rightEFunc,
        //    isLeftBorderEnabled, isRightBorderEnabled)
        //);
        generator.reset(new FDTD::FieldGeneratorType(this));
    }

    inline void FDTD::setTimeStep(FP dt)
    {
        if (ifCourantConditionSatisfied(dt)) {
            this->dt = dt;
            if (pml) pml.reset(new FDTD::PmlType(this, pml->sizePML));
            if (boundaryCondition) boundaryCondition.reset(boundaryCondition->createInstance(this));
            if (generator) generator.reset(new FDTD::FieldGeneratorType(this));  // TODO
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
        if (boundaryCondition) boundaryCondition->generateB();
        if (generator) generator->generateB();
        updateE();
        if (pml) pml->updateE();
        if (boundaryCondition) boundaryCondition->generateE();
        if (generator) generator->generateE();
        updateHalfB();
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

        const Int3 begin = internalBAreaBegin;
        const Int3 end = internalBAreaEnd;

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

        const Int3 begin = internalBAreaBegin;
        const Int3 end = internalBAreaEnd;

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

        const Int3 begin = internalBAreaBegin;
        const Int3 end = internalBAreaEnd;

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

        const Int3 begin = internalEAreaBegin;
        const Int3 end = internalEAreaEnd;

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

        const Int3 begin = internalEAreaBegin;
        const Int3 end = internalEAreaEnd;

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

        const Int3 begin = internalEAreaBegin;
        const Int3 end = internalEAreaEnd;

        OMP_FOR()
        for (int i = begin.x; i < end.x; i++) {
            grid->Ex(i, 0, 0) += coeffCurrent * grid->Jx(i, 0, 0);
            grid->Ey(i, 0, 0) += coeffCurrent * grid->Jy(i, 0, 0) -
                coeffXY * (grid->Bz(i + 1, 0, 0) - grid->Bz(i, 0, 0));
            grid->Ez(i, 0, 0) += coeffCurrent * grid->Jz(i, 0, 0) +
                coeffXZ * (grid->By(i + 1, 0, 0) - grid->By(i, 0, 0));
        }
    }

}
