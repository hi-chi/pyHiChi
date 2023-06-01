#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"
#include "PmlPstd.h"
#include "FieldBoundaryConditionSpectral.h"
#include "FieldGeneratorSpectral.h"

namespace pfc {
    class PSTD : public SpectralFieldSolver<PSTDGridType>
    {
    public:

        using GridType = PSTDGrid;
        using PmlType = PmlPstd;
        using FieldGeneratorType = FieldGeneratorSpectral<PSTDGridType>;
        using PeriodicalBoundaryConditionType = PeriodicalBoundaryConditionPstd;

        PSTD(PSTDGrid * grid, double dt);

        void updateFields();

        void updateHalfB();
        void updateE();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void setBoundaryCondition(
            FieldBoundaryCondition<GridTypes::PSTDGridType>* _boundaryCondition);

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
            double tmp = sqrt(1.0 / (grid->steps.x*grid->steps.x) +
                1.0 / (grid->steps.y*grid->steps.y) +
                1.0 / (grid->steps.z*grid->steps.z));
            return 2.0 / (constants::pi * constants::c * tmp);
        }

        bool ifCourantConditionSatisfied(FP dt) const {
            return dt < getCourantCondition();
        }

    private:

        PmlSpectral<GridTypes::PSTDGridType>* getPml() {
            return (PmlSpectral<GridTypes::PSTDGridType>*)pml.get();
        }

    };

    inline PSTD::PSTD(PSTDGrid* grid, double dt) :
        SpectralFieldSolver<GridTypes::PSTDGridType>(grid, dt)
    {
        if (!ifCourantConditionSatisfied(dt)) {
            std::cout
                << "WARNING: PSTD Courant condition is not satisfied. Another time step was setted up"
                << std::endl;
            this->dt = getCourantCondition() * 0.5;
        }
        updateDims();
        updateInternalDims();
    }

    inline void PSTD::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        pml.reset(new PmlPstd(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    inline void PSTD::setBoundaryCondition(
        FieldBoundaryCondition<GridTypes::PSTDGridType>* _boundaryCondition)
    {
        boundaryCondition.reset(_boundaryCondition->createInstance(this));
    }

    inline void PSTD::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        FieldGeneratorType::FunctionType bxFunc, FieldGeneratorType::FunctionType byFunc,
        FieldGeneratorType::FunctionType bzFunc, FieldGeneratorType::FunctionType exFunc,
        FieldGeneratorType::FunctionType eyFunc, FieldGeneratorType::FunctionType ezFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        //generator.reset(new PSTD::FieldGeneratorType(
        //    this, leftGenIndex, rightGenIndex,
        //    bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
        //    isLeftBorderEnabled, isRightBorderEnabled)
        //);
        generator.reset(new PSTD::FieldGeneratorType(this));
    }

    inline void PSTD::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftBFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightBFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftEFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightEFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        //generator.reset(new PSTD::FieldGeneratorType(
        //    this, leftGenIndex, rightGenIndex,
        //    leftBFunc, rightBFunc, leftEFunc, rightEFunc,
        //    isLeftBorderEnabled, isRightBorderEnabled)
        //);
        generator.reset(new PSTD::FieldGeneratorType(this));
    }

    inline void PSTD::setTimeStep(FP dt)
    {
        if (ifCourantConditionSatisfied(dt)) {
            this->dt = dt;
            if (pml) pml.reset(new PSTD::PmlType(this, pml->sizePML));
            if (boundaryCondition) boundaryCondition.reset(boundaryCondition->createInstance(this));
            if (generator) generator.reset(new PSTD::FieldGeneratorType(this));  // TODO
        }
        else {
            std::cout
                << "WARNING: PSTD Courant condition is not satisfied. Time step was not changed"
                << std::endl;
        }
    }

    inline void PSTD::updateFields()
    {
        doFourierTransform(fourier_transform::Direction::RtoC);

        if (pml) getPml()->updateBSplit();
        updateHalfB();

        if (pml) getPml()->updateESplit();
        updateE();

        if (pml) getPml()->updateBSplit();
        updateHalfB();

        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml) getPml()->doSecondStep();

        globalTime += dt;
    }

    inline void PSTD::updateHalfB()
    {
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;
        double dt = 0.5 * this->dt;
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
        const Int3 begin = updateComplexEAreaBegin;
        const Int3 end = updateComplexEAreaEnd;
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

}
