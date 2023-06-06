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
    class PSATDTimeStraggeredT : public SpectralFieldSolver<PSATDTimeStraggeredGridType>
    {
    public:

        using GridType = PSATDTimeStraggeredGrid;
        using PmlType = PmlPsatdTimeStraggered<PSATDTimeStraggeredGridType>;
        using FieldGeneratorType = FieldGeneratorSpectral<PSATDTimeStraggeredGridType>;
        using PeriodicalBoundaryConditionType = PeriodicalBoundaryConditionSpectral<PSATDTimeStraggeredGridType>;

        PSATDTimeStraggeredT(GridType* grid, FP dt);

        void updateFields();

        void updateHalfB();
        void updateE();

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

        ScalarField<complexFP> tmpJx, tmpJy, tmpJz;

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

        void save(std::ostream& ostr);
        void load(std::istream& istr);

    protected:

        PmlType* getPml() const {
            return static_cast<PmlType*>(pml.get());
        }
        FieldGeneratorType* getGenerator() const {
            return static_cast<FieldGeneratorType*>(generator.get());
        }

        void saveJ();
        void assignJ(SpectralScalarField<FP, complexFP>& J, ScalarField<complexFP>& tmpJ);
    };

    template<bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::save(std::ostream& ostr)
    {
        FieldSolver<GridType::gridType>::save(ostr);
        tmpJx.save(ostr);
        tmpJy.save(ostr);
        tmpJz.save(ostr);
    }

    template<bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::load(std::istream& istr)
    {
        FieldSolver<GridType::gridType>::load(istr);
        tmpJx.load(istr);
        tmpJy.load(istr);
        tmpJz.load(istr);
    }

    template <bool ifPoisson>
    inline PSATDTimeStraggeredT<ifPoisson>::PSATDTimeStraggeredT(GridType* _grid, FP dt) :
        SpectralFieldSolver<GridType::gridType>(_grid, dt),
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
        pml.reset(new PmlType(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        updateInternalDims();
    }

    template <bool ifPoisson>
    template <class TBoundaryCondition>
    inline void PSATDTimeStraggeredT<ifPoisson>::setBoundaryCondition()
    {
        for (int d = 0; d < grid->dimensionality; d++)
            boundaryConditions[d].reset(new TBoundaryCondition((CoordinateEnum)d, this));
    }

    template <bool ifPoisson>
    template <class TBoundaryCondition>
    inline void PSATDTimeStraggeredT<ifPoisson>::setBoundaryCondition(CoordinateEnum axis)
    {
        if ((int)axis >= grid->dimensionality) {
            std::cout
                << "WARNING: an attempt to set boundary conditions for an axis greater than the dimensionality is ignored"
                << std::endl;
        }
        boundaryConditions[(int)axis].reset(new TBoundaryCondition(axis, this));
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        FieldGeneratorType::FunctionType bxFunc, FieldGeneratorType::FunctionType byFunc,
        FieldGeneratorType::FunctionType bzFunc, FieldGeneratorType::FunctionType exFunc,
        FieldGeneratorType::FunctionType eyFunc, FieldGeneratorType::FunctionType ezFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        generator.reset(new FieldGeneratorType(
            this, leftGenIndex, rightGenIndex,
            bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
        );
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftBFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightBFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& leftEFunc,
        const std::array<std::array<FieldGeneratorType::FunctionType, 3>, 3>& rightEFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        generator.reset(new FieldGeneratorType(
            this, leftGenIndex, rightGenIndex,
            leftBFunc, rightBFunc, leftEFunc, rightEFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
        );
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::setTimeStep(FP dt)
    {
        this->dt = dt;
        if (pml) pml.reset(new PmlType(this, pml->sizePML));
        for (int d = 0; d < 3; d++)
            if (boundaryConditions[d])
                boundaryConditions[d].reset(boundaryConditions[d]->createInstance(this));
        if (generator) generator.reset(new FieldGeneratorType(*getGenerator()));
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::assignJ(
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

        if (pml) getPml()->updateBSplit();
        updateHalfB();
        for (int d = 0; d < grid->dimensionality; d++)
            if (boundaryConditions[d]) boundaryConditions[d]->generateB();

        if (pml) getPml()->updateESplit();
        updateE();
        for (int d = 0; d < grid->dimensionality; d++)
            if (boundaryConditions[d]) boundaryConditions[d]->generateE();

        if (pml) getPml()->updateBSplit();
        updateHalfB();
        for (int d = 0; d < grid->dimensionality; d++)
            if (boundaryConditions[d]) boundaryConditions[d]->generateB();

        saveJ();
        doFourierTransform(fourier_transform::Direction::CtoR);

        if (pml) getPml()->updateB();
        if (pml) getPml()->updateE();

        globalTime += dt;
    }

    template <bool ifPoisson>
    inline void PSATDTimeStraggeredT<ifPoisson>::convertFieldsPoissonEquation() {
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
    inline void PSATDTimeStraggeredT<ifPoisson>::updateHalfB()
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
    inline void PSATDTimeStraggeredT<ifPoisson>::updateE()
    {
        const Int3 begin = updateComplexEAreaBegin;
        const Int3 end = updateComplexEAreaEnd;

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
    inline void PSATDTimeStraggeredT<true>::updateE()
    {
        const Int3 begin = updateComplexEAreaBegin;
        const Int3 end = updateComplexEAreaEnd;

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

    typedef PSATDTimeStraggeredT<true> PSATDTimeStraggeredPoisson;
    typedef PSATDTimeStraggeredT<false> PSATDTimeStraggered;

}
