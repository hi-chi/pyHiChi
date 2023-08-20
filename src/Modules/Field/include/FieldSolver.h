#pragma once
#include "Grid.h"
#include "SpectralGrid.h"
#include "FourierTransform.h"

#include <algorithm>
#include <memory>
#include <array>
#include <typeinfo>

namespace pfc {
    template<class GridType>
    class FieldBoundaryCondition;

    // Base class for field solvers on Grid
    template<class SchemeParams>
    class FieldSolver
    {
    public:

        FieldSolver(typename SchemeParams::GridType* grid, FP dt);

        // constructor for loading
        explicit FieldSolver(typename SchemeParams::GridType* grid);

        /* implement the next methods in derived classes
        void updateFields();
        void setTimeStep(FP dt);
        static FP getCourantConditionTimeStep(const FP3& gridSteps);
        FP getCourantConditionTimeStep() const;
        bool isCourantConditionSatisfied(FP dt) const;
        static bool isTimeStaggered();
        void save(std::ostream& ostr);
        void load(std::istream& istr);
        */

        FP getTimeStep() const { return this->dt; }

        void setTime(FP t) { this->globalTime = t; }
        FP getTime() const { return this->globalTime; }

        void updateDomainBorders();

        void applyBoundaryConditionsB(FP time);
        void applyBoundaryConditionsE(FP time);

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        void setFieldGenerator(
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            typename SchemeParams::FieldGeneratorType::FunctionType bxFunc, typename SchemeParams::FieldGeneratorType::FunctionType byFunc,
            typename SchemeParams::FieldGeneratorType::FunctionType bzFunc, typename SchemeParams::FieldGeneratorType::FunctionType exFunc,
            typename SchemeParams::FieldGeneratorType::FunctionType eyFunc, typename SchemeParams::FieldGeneratorType::FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));
        void setFieldGenerator(
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& xLeftBFunc,  // { bx, by, bz }
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& xRightBFunc,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& yLeftBFunc,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& yRightBFunc,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& zLeftBFunc,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& zRightBFunc,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& xLeftEFunc,  // { ex, ey, ez }
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& xRightEFunc,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& yLeftEFunc,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& yRightEFunc,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& zLeftEFunc,
            const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& zRightEFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1));
        void saveFieldGenerator(std::ostream& ostr);
        void loadFieldGenerator(std::istream& istr);
        void resetFieldGenerator();

        /* next methods are already implemented in derived classes
        void setPML(Int3 sizePML);
        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void savePML(std::ostream& ostr);
        void loadPML(std::istream& istr);
        void resetPML();
        */

        /* implement the next methods in derived classes
        void set<Type>BoundaryConditions();
        void set<Type>BoundaryConditions(CoordinateEnum axis);
        void saveBoundaryConditions(std::ostream& ostr);
        void loadBoundaryConditions(std::istream& istr);
        */
        void resetBoundaryConditions();

        typename SchemeParams::GridType* grid = nullptr;

        std::unique_ptr<typename SchemeParams::PmlType> pml;
        std::array<std::unique_ptr<FieldBoundaryCondition<typename SchemeParams::GridType>>, 3> boundaryConditions;
        std::unique_ptr<typename SchemeParams::FieldGeneratorType> generator;

        // Index space being updated in non-pml area
        Int3 internalIndexBegin, internalIndexEnd;
        // Internal index space, const
        const Int3 domainIndexBegin, domainIndexEnd;

        FP globalTime = 0.0;
        FP dt = 0.0;
    };

    template<class SchemeParams>
    inline FieldSolver<SchemeParams>::FieldSolver(typename SchemeParams::GridType* grid, FP dt) :
        grid(grid), globalTime(0.0), dt(dt), domainIndexBegin(grid->getNumExternalLeftCells()),
        domainIndexEnd(grid->getNumExternalLeftCells() + grid->numInternalCells)
    {
        this->updateDomainBorders();
    }

    template<class SchemeParams>
    inline FieldSolver<SchemeParams>::FieldSolver(typename SchemeParams::GridType* grid) :
        grid(grid), domainIndexBegin(grid->getNumExternalLeftCells()),
        domainIndexEnd(grid->getNumExternalLeftCells() + grid->numInternalCells)
    {
        this->updateDomainBorders();
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::save(std::ostream& ostr)
    {
        ostr.write((char*)&globalTime, sizeof(globalTime));
        ostr.write((char*)&dt, sizeof(dt));
        // implement module saving in derived classes
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::load(std::istream& istr)
    {
        istr.read((char*)&globalTime, sizeof(globalTime));
        istr.read((char*)&dt, sizeof(dt));
        // implement module loading in derived classes
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::updateDomainBorders()
    {
        this->internalIndexBegin = this->grid->getNumExternalLeftCells();
        this->internalIndexEnd = this->internalIndexBegin + this->grid->numInternalCells;

        if (pml)
        {
            this->internalIndexBegin += pml->sizePML;
            this->internalIndexEnd -= pml->sizePML;
        }
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::applyBoundaryConditionsB(FP time)
    {
        for (int d = 0; d < grid->dimensionality; d++)
            if (boundaryConditions[d])
                boundaryConditions[d]->generateB(time);
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::applyBoundaryConditionsE(FP time)
    {
        for (int d = 0; d < grid->dimensionality; d++)
            if (boundaryConditions[d])
                boundaryConditions[d]->generateE(time);
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        typename SchemeParams::FieldGeneratorType::FunctionType bxFunc,
        typename SchemeParams::FieldGeneratorType::FunctionType byFunc,
        typename SchemeParams::FieldGeneratorType::FunctionType bzFunc,
        typename SchemeParams::FieldGeneratorType::FunctionType exFunc,
        typename SchemeParams::FieldGeneratorType::FunctionType eyFunc,
        typename SchemeParams::FieldGeneratorType::FunctionType ezFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        this->generator.reset(new typename SchemeParams::FieldGeneratorType(
            this->grid, this->dt, this->domainIndexBegin, this->domainIndexEnd,
            leftGenIndex, rightGenIndex,
            bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
        );
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::setFieldGenerator(
        const Int3& leftGenIndex, const Int3& rightGenIndex,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& xLeftBFunc,  // { bx, by, bz }
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& xRightBFunc,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& yLeftBFunc,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& yRightBFunc,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& zLeftBFunc,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& zRightBFunc,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& xLeftEFunc,  // { ex, ey, ez }
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& xRightEFunc,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& yLeftEFunc,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& yRightEFunc,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& zLeftEFunc,
        const std::array<typename SchemeParams::FieldGeneratorType::FunctionType, 3>& zRightEFunc,
        const Int3& isLeftBorderEnabled, const Int3& isRightBorderEnabled)
    {
        this->generator.reset(new typename SchemeParams::FieldGeneratorType(
            this->grid, this->dt, this->domainIndexBegin, this->domainIndexEnd,
            leftGenIndex, rightGenIndex,
            xLeftBFunc, xRightBFunc, yLeftBFunc, yRightBFunc, zLeftBFunc, zRightBFunc,
            xLeftEFunc, xRightEFunc, yLeftEFunc, yRightEFunc, zLeftEFunc, zRightEFunc,
            isLeftBorderEnabled, isRightBorderEnabled)
        );
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::saveFieldGenerator(std::ostream& ostr)
    {
        int isGenerator = this->generator ? 1 : 0;
        ostr.write((char*)&isGenerator, sizeof(isGenerator));
        if (isGenerator) this->generator->save(ostr);
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::loadFieldGenerator(std::istream& istr)
    {
        int isGenerator = 0;
        istr.read((char*)&isGenerator, sizeof(isGenerator));
        if (isGenerator) {
            this->generator.reset(new typename SchemeParams::FieldGeneratorType(
                this->grid, this->dt, this->domainIndexBegin, this->domainIndexEnd));
            this->generator->load(istr);
        }
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::resetFieldGenerator()
    {
        if (this->generator) this->generator.reset(new typename SchemeParams::FieldGeneratorType(this->grid, this->dt,
            this->domainIndexBegin, this->domainIndexEnd, *(this->generator)));
    }

    template<class SchemeParams>
    inline void FieldSolver<SchemeParams>::resetBoundaryConditions()
    {
        for (int d = 0; d < 3; d++)
            if (this->boundaryConditions[d])
                this->boundaryConditions[d].reset(this->boundaryConditions[d]->createInstance(this->grid,
                    this->domainIndexBegin, this->domainIndexEnd, this->boundaryConditions[d]->axis));
    }


    template<class SchemeParams>
    class RealFieldSolver : public FieldSolver<SchemeParams>
    {
    public:

        RealFieldSolver(typename SchemeParams::GridType* grid, FP dt) :
            FieldSolver<SchemeParams>(grid, dt)
        {}

        // constructor for loading
        explicit RealFieldSolver(typename SchemeParams::GridType* grid) :
            FieldSolver<SchemeParams>(grid)
        {}

        void setPML(Int3 sizePML);
        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void savePML(std::ostream& ostr);
        void loadPML(std::istream& istr);
        void resetPML();

    private:
        // Copy and assignment are disallowed.
        RealFieldSolver(const RealFieldSolver&);
        RealFieldSolver& operator =(const RealFieldSolver&);
    };

    template<class SchemeParams>
    inline void RealFieldSolver<SchemeParams>::setPML(Int3 sizePML)
    {
        this->pml.reset(new typename SchemeParams::PmlType(this->grid, this->dt,
            this->domainIndexBegin, this->domainIndexEnd, sizePML));
        this->updateDomainBorders();
    }

    template<class SchemeParams>
    inline void RealFieldSolver<SchemeParams>::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        this->setPML(Int3(sizePMLx, sizePMLy, sizePMLz));
    }

    template<class SchemeParams>
    inline void RealFieldSolver<SchemeParams>::savePML(std::ostream& ostr)
    {
        int isPml = this->pml ? 1 : 0;
        ostr.write((char*)&isPml, sizeof(isPml));
        if (isPml) this->pml->save(ostr);
    }

    template<class SchemeParams>
    inline void RealFieldSolver<SchemeParams>::loadPML(std::istream& istr)
    {
        int isPml = 0;
        istr.read((char*)&isPml, sizeof(isPml));
        if (isPml) {
            this->pml.reset(new typename SchemeParams::PmlType(this->grid, this->dt,
                this->domainIndexBegin, this->domainIndexEnd));
            this->pml->load(istr);
        }
    }

    template<class SchemeParams>
    inline void RealFieldSolver<SchemeParams>::resetPML()
    {
        if (this->pml) this->setPML(this->pml->sizePML);
    }


    template<class SchemeParams>
    class SpectralFieldSolver : public FieldSolver<SchemeParams>
    {
    public:
        SpectralFieldSolver(typename SchemeParams::GridType* grid, FP dt);

        // constructor for loading
        explicit SpectralFieldSolver(typename SchemeParams::GridType* grid);

        void doFourierTransformB(fourier_transform::Direction direction);
        void doFourierTransformE(fourier_transform::Direction direction);
        void doFourierTransformJ(fourier_transform::Direction direction);
        void doFourierTransform(fourier_transform::Direction direction);

        FP3 getWaveVector(const Int3& ind);

        void initComplexPart();
        void updateComplexDomainBorders();

        void setPML(Int3 sizePML);
        void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
        void savePML(std::ostream& ostr);
        void loadPML(std::istream& istr);
        void resetPML();

        std::unique_ptr<SpectralGrid<FP, complexFP>> complexGrid;

        Int3 complexDomainIndexBegin, complexDomainIndexEnd;

        FourierTransformGrid fourierTransform;

    private:
        // Copy and assignment are disallowed.
        SpectralFieldSolver(const SpectralFieldSolver&);
        SpectralFieldSolver& operator=(const SpectralFieldSolver&);
    };

    template<class SchemeParams>
    inline SpectralFieldSolver<SchemeParams>::SpectralFieldSolver(typename SchemeParams::GridType* grid, FP dt) :
        FieldSolver<SchemeParams>(grid, dt)
    {
        initComplexPart();
    }

    template<class SchemeParams>
    inline SpectralFieldSolver<SchemeParams>::SpectralFieldSolver(typename SchemeParams::GridType* grid) :
        FieldSolver<SchemeParams>(grid)
    {
        initComplexPart();
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::doFourierTransformB(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::x, direction);
        fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::y, direction);
        fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::z, direction);
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::doFourierTransformE(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(FieldEnum::E, CoordinateEnum::x, direction);
        fourierTransform.doFourierTransform(FieldEnum::E, CoordinateEnum::y, direction);
        fourierTransform.doFourierTransform(FieldEnum::E, CoordinateEnum::z, direction);
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::doFourierTransformJ(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(FieldEnum::J, CoordinateEnum::x, direction);
        fourierTransform.doFourierTransform(FieldEnum::J, CoordinateEnum::y, direction);
        fourierTransform.doFourierTransform(FieldEnum::J, CoordinateEnum::z, direction);
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::doFourierTransform(fourier_transform::Direction direction)
    {
        doFourierTransformE(direction);
        doFourierTransformB(direction);
        doFourierTransformJ(direction);
    }

    template<class SchemeParams>
    inline FP3 SpectralFieldSolver<SchemeParams>::getWaveVector(const Int3& ind)
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

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::initComplexPart() {
        this->complexGrid.reset(new SpectralGrid<FP, complexFP>(this->grid,
            fourier_transform::getSizeOfComplexArray(this->domainIndexEnd - this->domainIndexBegin)
            ));
        this->fourierTransform.template initialize<typename SchemeParams::GridType>(this->grid, this->complexGrid.get());
        this->updateComplexDomainBorders();
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::updateComplexDomainBorders()
    {
        this->complexDomainIndexBegin = Int3(0, 0, 0);
        this->complexDomainIndexEnd = this->complexGrid->numCells;
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::setPML(Int3 sizePML)
    {
        this->pml.reset(new typename SchemeParams::PmlType(this->grid, this->complexGrid.get(), this->dt,
            this->domainIndexBegin, this->domainIndexEnd,
            this->complexDomainIndexBegin, this->complexDomainIndexEnd, sizePML));
        this->updateDomainBorders();
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
    {
        this->setPML(Int3(sizePMLx, sizePMLy, sizePMLz));
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::savePML(std::ostream& ostr)
    {
        int isPml = this->pml ? 1 : 0;
        ostr.write((char*)&isPml, sizeof(isPml));
        if (isPml) this->pml->save(ostr);
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::loadPML(std::istream& istr)
    {
        int isPml = 0;
        istr.read((char*)&isPml, sizeof(isPml));
        if (isPml) {
            this->pml.reset(new typename SchemeParams::PmlType(this->grid, this->complexGrid.get(), this->dt,
                this->domainIndexBegin, this->domainIndexEnd,
                this->complexDomainIndexBegin, this->complexDomainIndexEnd));
            this->pml->load(istr);
        }
    }

    template<class SchemeParams>
    inline void SpectralFieldSolver<SchemeParams>::resetPML()
    {
        if (this->pml) this->setPML(this->pml->sizePML);
    }

};
