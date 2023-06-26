#pragma once
#include "Grid.h"
#include "SpectralGrid.h"
#include "FourierTransform.h"

#include <algorithm>
#include <memory>
#include <array>

namespace pfc {
    template<class TGrid>
    class FieldBoundaryCondition;

    // Base class for field solvers on Grid
    template<class TGrid, class TPml, class TFieldGenerator>
    class FieldSolver
    {

    public:

        FieldSolver(TGrid* grid);  // use when load
        FieldSolver(TGrid* grid, FP dt);

        /* implement the next methods in derived classes
        void updateFields();
        void setTimeStep(FP dt);
        static FP getCourantConditionTimeStep(const FP3& gridSteps);
        FP getCourantConditionTimeStep() const;
        bool isCourantConditionSatisfied(FP dt) const;
        static bool isTimeStaggered();
        void save(std::ostream& ostr);  // of necessity
        void load(std::istream& istr);  // of necessity
        */

        FP getTimeStep() const { return this->dt; }

        void setTime(FP t) { this->globalTime = t; }
        FP getTime() const { return this->globalTime; }

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        void updateDomainBorders();

        void applyBoundaryConditionsB(FP time);
        void applyBoundaryConditionsE(FP time);

        void setFieldGenerator(
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            typename TFieldGenerator::FunctionType bxFunc, typename TFieldGenerator::FunctionType byFunc,
            typename TFieldGenerator::FunctionType bzFunc, typename TFieldGenerator::FunctionType exFunc,
            typename TFieldGenerator::FunctionType eyFunc, typename TFieldGenerator::FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));
        void setFieldGenerator(
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            /* first index is index of edge (x, y, z),
            second index is index of field component (ex, ey, ez or bx, by, bz) */
            const std::array<std::array<typename TFieldGenerator::FunctionType, 3>, 3>& leftBFunc,
            const std::array<std::array<typename TFieldGenerator::FunctionType, 3>, 3>& rightBFunc,
            const std::array<std::array<typename TFieldGenerator::FunctionType, 3>, 3>& leftEFunc,
            const std::array<std::array<typename TFieldGenerator::FunctionType, 3>, 3>& rightEFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));
        void resetFieldGenerator();

        /* next methods are implemented in derived classes
        void setPML(Int3 sizePML);
        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void resetPML();
        */

        template <class TBoundaryCondition>
        void setBoundaryCondition();
        template <class TBoundaryCondition>
        void setBoundaryCondition(CoordinateEnum axis);
        void resetBoundaryConditions();

        TGrid* grid = nullptr;

        std::unique_ptr<TPml> pml;
        std::array<std::unique_ptr<FieldBoundaryCondition<TGrid>>, 3> boundaryConditions;
        std::unique_ptr<TFieldGenerator> generator;

        // Index space being updated in non-pml area
        Int3 internalIndexBegin, internalIndexEnd;
        // Internal index space, const
        const Int3 domainIndexBegin, domainIndexEnd;

        FP globalTime = 0.0;
        FP dt = 0.0;
    };

    template<class TGrid, class TPml, class TFieldGenerator>
    inline FieldSolver<TGrid, TPml, TFieldGenerator>::FieldSolver(TGrid* grid) :
        grid(grid), domainIndexBegin(grid->getNumExternalLeftCells()),
        domainIndexEnd(grid->getNumExternalLeftCells() + grid->numInternalCells)
    {
        updateDomainBorders();
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline FieldSolver<TGrid, TPml, TFieldGenerator>::FieldSolver(TGrid* grid, FP dt) :
        grid(grid), globalTime(0.0), dt(dt), domainIndexBegin(grid->getNumExternalLeftCells()),
        domainIndexEnd(grid->getNumExternalLeftCells() + grid->numInternalCells)
    {
        updateDomainBorders();
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::save(std::ostream& ostr)
    {
        ostr.write((char*)&globalTime, sizeof(globalTime));
        ostr.write((char*)&dt, sizeof(dt));
        // TODO: save modules
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::load(std::istream& istr)
    {
        istr.read((char*)&globalTime, sizeof(globalTime));
        istr.read((char*)&dt, sizeof(dt));
        // TODO: load modules
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::updateDomainBorders()
    {
        this->internalIndexBegin = this->grid->getNumExternalLeftCells();
        this->internalIndexEnd = this->internalIndexBegin + this->grid->numInternalCells;

        if (pml)
        {
            this->internalIndexBegin += pml->sizePML;
            this->internalIndexEnd -= pml->sizePML;
        }
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::applyBoundaryConditionsB(FP time)
    {
        for (int d = 0; d < grid->dimensionality; d++)
            if (boundaryConditions[d])
                boundaryConditions[d]->generateB(time);
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::applyBoundaryConditionsE(FP time)
    {
        for (int d = 0; d < grid->dimensionality; d++)
            if (boundaryConditions[d])
                boundaryConditions[d]->generateE(time);
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        typename TFieldGenerator::FunctionType bxFunc, typename TFieldGenerator::FunctionType byFunc,
        typename TFieldGenerator::FunctionType bzFunc, typename TFieldGenerator::FunctionType exFunc,
        typename TFieldGenerator::FunctionType eyFunc, typename TFieldGenerator::FunctionType ezFunc,
        const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
        const Int3& isRightBorderEnabled = Int3(1, 1, 1))
    {
        generator.reset(new TFieldGenerator(
            this->grid, this->dt, leftGenIndex, rightGenIndex,
            bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
        );
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        /* first index is index of edge (x, y, z),
        second index is index of field component (ex, ey, ez or bx, by, bz) */
        const std::array<std::array<typename TFieldGenerator::FunctionType, 3>, 3>& leftBFunc,
        const std::array<std::array<typename TFieldGenerator::FunctionType, 3>, 3>& rightBFunc,
        const std::array<std::array<typename TFieldGenerator::FunctionType, 3>, 3>& leftEFunc,
        const std::array<std::array<typename TFieldGenerator::FunctionType, 3>, 3>& rightEFunc,
        const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
        const Int3& isRightBorderEnabled = Int3(1, 1, 1))
    {
        generator.reset(new TFieldGenerator(
            this->grid, this->dt, leftGenIndex, rightGenIndex,
            leftBFunc, rightBFunc, leftEFunc, rightEFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
        );
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::resetFieldGenerator()
    {
        if (generator) generator.reset(new TFieldGenerator(this->grid, this->dt, *generator));
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    template <class TBoundaryCondition>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::setBoundaryCondition()
    {
        for (int d = 0; d < this->grid->dimensionality; d++)
            boundaryConditions[d].reset(new TBoundaryCondition(this->grid,
                (CoordinateEnum)d, this->domainIndexBegin, this->domainIndexEnd));
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    template <class TBoundaryCondition>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::setBoundaryCondition(CoordinateEnum axis)
    {
        if ((int)axis >= this->grid->dimensionality) {
            std::cout
                << "WARNING: an attempt to set boundary conditions for an axis greater than the dimensionality is ignored"
                << std::endl;
        }
        boundaryConditions[(int)axis].reset(new TBoundaryCondition(this->grid,
            axis, this->domainIndexBegin, this->domainIndexEnd));
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void FieldSolver<TGrid, TPml, TFieldGenerator>::resetBoundaryConditions()
    {
        for (int d = 0; d < 3; d++)
            if (boundaryConditions[d])
                boundaryConditions[d].reset(boundaryConditions[d]->createInstance(this->grid,
                    boundaryConditions[d]->axis, this->domainIndexBegin, this->domainIndexEnd));
    }


    template<class TGrid, class TPml, class TFieldGenerator>
    class RealFieldSolver : public FieldSolver<TGrid, TPml, TFieldGenerator>
    {
    public:
        RealFieldSolver(TGrid* grid) :  // use when load
            FieldSolver<TGrid, TPml, TFieldGenerator>(grid)
        {}
        RealFieldSolver(TGrid* _grid, FP dt) :
            FieldSolver<TGrid, TPml, TFieldGenerator>(_grid, dt)
        {}

        void setPML(Int3 sizePML);
        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void resetPML();

    private:
        // Copy and assignment are disallowed.
        RealFieldSolver(const RealFieldSolver&);
        RealFieldSolver& operator =(const RealFieldSolver&);
    };

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void RealFieldSolver<TGrid, TPml, TFieldGenerator>::setPML(Int3 sizePML)
    {
        pml.reset(new TPml(this->grid, this->dt, sizePML,
            this->domainIndexBegin, this->domainIndexEnd));
        updateDomainBorders();
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void RealFieldSolver<TGrid, TPml, TFieldGenerator>::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        setPML(Int3(sizePMLx, sizePMLy, sizePMLz));
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void RealFieldSolver<TGrid, TPml, TFieldGenerator>::resetPML()
    {
        if (pml) setPML(pml->sizePML);
    }


    template<class TGrid, class TPml, class TFieldGenerator>
    class SpectralFieldSolver : public FieldSolver<TGrid, TPml, TFieldGenerator>
    {
    public:
        SpectralFieldSolver(TGrid* grid);  // use when load
        SpectralFieldSolver(TGrid* grid, FP dt);

        void doFourierTransformB(fourier_transform::Direction direction);
        void doFourierTransformE(fourier_transform::Direction direction);
        void doFourierTransformJ(fourier_transform::Direction direction);
        void doFourierTransform(fourier_transform::Direction direction);

        FP3 getWaveVector(const Int3& ind);

        void initComplexPart();
        void updateComplexDomainBorders();

        void setPML(Int3 sizePML);
        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void resetPML();

        std::unique_ptr<SpectralGrid<FP, complexFP>> complexGrid;

        Int3 complexDomainIndexBegin, complexDomainIndexEnd;

        FourierTransformGrid fourierTransform;

    private:
        // Copy and assignment are disallowed.
        SpectralFieldSolver(const SpectralFieldSolver&);
        SpectralFieldSolver& operator=(const SpectralFieldSolver&);
    };

    template<class TGrid, class TPml, class TFieldGenerator>
    inline SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::SpectralFieldSolver(TGrid* grid) :
        FieldSolver<TGrid, TPml, TFieldGenerator>(grid)
    {
        initComplexPart();
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::SpectralFieldSolver(TGrid* grid, FP dt) :
        FieldSolver<TGrid, TPml, TFieldGenerator>(grid, dt)
    {
        initComplexPart();
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::doFourierTransformB(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::x, direction);
        fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::y, direction);
        fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::z, direction);
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::doFourierTransformE(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(FieldEnum::E, CoordinateEnum::x, direction);
        fourierTransform.doFourierTransform(FieldEnum::E, CoordinateEnum::y, direction);
        fourierTransform.doFourierTransform(FieldEnum::E, CoordinateEnum::z, direction);
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::doFourierTransformJ(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(FieldEnum::J, CoordinateEnum::x, direction);
        fourierTransform.doFourierTransform(FieldEnum::J, CoordinateEnum::y, direction);
        fourierTransform.doFourierTransform(FieldEnum::J, CoordinateEnum::z, direction);
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::doFourierTransform(fourier_transform::Direction direction)
    {
        doFourierTransformE(direction);
        doFourierTransformB(direction);
        doFourierTransformJ(direction);
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline FP3 SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::getWaveVector(const Int3& ind)
    {
        const Int3 domainSize = this->domainIndexEnd - this->domainIndexBegin;
        FP kx = ((FP)2.0 * constants::pi * ((ind.x <= domainSize.x / 2) ? ind.x : ind.x - domainSize.x)) /
            (this->grid->steps.x * domainSize.x);
        FP ky = ((FP)2.0 * constants::pi * ((ind.y <= domainSize.y / 2) ? ind.y : ind.y - domainSize.y)) /
            (this->grid->steps.y * domainSize.y);
        FP kz = ((FP)2.0 * constants::pi * ((ind.z <= domainSize.z / 2) ? ind.z : ind.z - domainSize.z)) /
            (this->grid->steps.z * domainSize.z);
        return FP3(kx, ky, kz);
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::initComplexPart() {
        complexGrid.reset(new SpectralGrid<FP, complexFP>(
            fourier_transform::getSizeOfComplexArray(this->grid->numCells),
            this->grid->globalGridDims, this->grid));
        fourierTransform.initialize<TGrid>(this->grid, this->complexGrid.get());
        updateComplexDomainBorders();
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::updateComplexDomainBorders()
    {
        this->complexDomainIndexBegin = Int3(0, 0, 0);
        this->complexDomainIndexEnd = this->complexGrid->numCells;
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::setPML(Int3 sizePML)
    {
        pml.reset(new TPml(this->grid, this->complexGrid.get(), this->dt, sizePML,
            this->domainIndexBegin, this->domainIndexEnd,
            this->complexDomainIndexBegin, this->complexDomainIndexEnd));
        updateDomainBorders();
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        setPML(Int3(sizePMLx, sizePMLy, sizePMLz));
    }

    template<class TGrid, class TPml, class TFieldGenerator>
    inline void SpectralFieldSolver<TGrid, TPml, TFieldGenerator>::resetPML()
    {
        if (pml) setPML(pml->sizePML);
    }

};
