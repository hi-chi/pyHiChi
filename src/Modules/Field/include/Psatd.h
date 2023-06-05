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
    class PSATDT : public SpectralFieldSolver<PSATDGridType>
    {
    public:

        using GridType = PSATDGrid;
        using PmlType = PmlPsatd;
        using FieldGeneratorType = FieldGeneratorSpectral<PSATDGridType>;
        using PeriodicalBoundaryConditionType = PeriodicalBoundaryConditionPsatd;

        PSATDT(GridType* grid, FP dt);

        void updateFields();

        virtual void updateEB();

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);

        template <class TBoundaryCondition>
        void setBoundaryCondition();
        template <class TBoundaryCondition>
        void setBoundaryCondition(CoordinateEnum axis);

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

        void convertFieldsPoissonEquation();

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

        bool ifCourantConditionSatisfied(FP dt) const {
            return true;
        }

    private:

        PmlType* getPml() const {
            return static_cast<PmlType*>(pml.get());
        }
        FieldGeneratorType* getGenerator() const {
            return static_cast<FieldGeneratorType*>(generator.get());
        }

    };

    template <bool ifPoisson>
    inline PSATDT<ifPoisson>::PSATDT(PSATDT<ifPoisson>::GridType* _grid, FP dt) :
        SpectralFieldSolver<GridTypes::PSATDGridType>(_grid, dt)
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
    template <class TBoundaryCondition>
    inline void PSATDT<ifPoisson>::setBoundaryCondition()
    {
        for (int d = 0; d < grid->dimensionality; d++)
            boundaryConditions[d].reset(new TBoundaryCondition((CoordinateEnum)d, this));
    }

    template <bool ifPoisson>
    template <class TBoundaryCondition>
    inline void PSATDT<ifPoisson>::setBoundaryCondition(CoordinateEnum axis)
    {
        if ((int)axis >= grid->dimensionality) {
            std::cout
                << "WARNING: an attempt to set boundary conditions for an axis greater than the dimensionality is ignored"
                << std::endl;
        }
        boundaryConditions[(int)axis].reset(new TBoundaryCondition(axis, this));
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        FieldGeneratorType::FunctionType bxFunc, FieldGeneratorType::FunctionType byFunc,
        FieldGeneratorType::FunctionType bzFunc, FieldGeneratorType::FunctionType exFunc,
        FieldGeneratorType::FunctionType eyFunc, FieldGeneratorType::FunctionType ezFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        generator.reset(new PSATDT<ifPoisson>::FieldGeneratorType(
            this, leftGenIndex, rightGenIndex,
            bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
        );
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftBFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightBFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftEFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightEFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        generator.reset(new PSATDT<ifPoisson>::::FieldGeneratorType(
            this, leftGenIndex, rightGenIndex,
            leftBFunc, rightBFunc, leftEFunc, rightEFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
        );
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::setTimeStep(FP dt)
    {
        this->dt = dt;
        if (pml) pml.reset(new PSATDT<ifPoisson>::PmlType(this, pml->sizePML));
        for (int d = 0; d < 3; d++)
            if (boundaryConditions[d])
                boundaryConditions[d].reset(boundaryConditions[d]->createInstance(this));
        if (generator) generator.reset(new PSATDT<ifPoisson>::FieldGeneratorType(*getGenerator()));
    }

    template <bool ifPoisson>
    inline void PSATDT<ifPoisson>::updateFields() {
        doFourierTransform(fourier_transform::Direction::RtoC);

        if (pml) getPml()->updateBSplit();
        updateEB();
        for (int d = 0; d < grid->dimensionality; d++)
            if (boundaryConditions[d]) {
                boundaryConditions[d]->generateE();
                boundaryConditions[d]->generateB();
            }

        if (pml) getPml()->updateESplit();
        updateEB();
        for (int d = 0; d < grid->dimensionality; d++)
            if (boundaryConditions[d]) {
                boundaryConditions[d]->generateE();
                boundaryConditions[d]->generateB();
            }

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
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;

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
        const Int3 begin = updateComplexBAreaBegin;
        const Int3 end = updateComplexBAreaEnd;

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

    typedef PSATDT<true> PSATDPoisson;
    typedef PSATDT<false> PSATD;

}
