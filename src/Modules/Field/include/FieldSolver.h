#pragma once
#include "Grid.h"
#include "SpectralGrid.h"
#include "Vectors.h"
#include "FourierTransform.h"

#include <algorithm>
#include <memory>
#include <array>

namespace pfc {
    template<GridTypes gridType>
    class Pml;
    template<GridTypes gridType>
    class FieldBoundaryCondition;
    template<GridTypes gridType>
    class FieldGenerator;

    // Base class for field solvers on Grid
    template<GridTypes gridType>
    class FieldSolver
    {

    public:

        FieldSolver(Grid<FP, gridType>* _grid, FP dt) : dt(dt)
        {
            this->grid = _grid;
            this->globalTime = (FP)0.0;
            this->generator.reset(nullptr);
        }

        virtual ~FieldSolver() {}
        
        void updateDims();
        void updateInternalDims();

        virtual void save(std::ostream& ostr);
        virtual void load(std::istream& istr);

        Grid<FP, gridType>* grid;

        std::unique_ptr<Pml<gridType>> pml;
        std::array<std::unique_ptr<FieldBoundaryCondition<gridType>>, 3> boundaryConditions;
        std::unique_ptr<FieldGenerator<gridType>> generator;

        // Index space being updated in form [begin, end).
        Int3 updateBAreaBegin, updateBAreaEnd;
        Int3 updateEAreaBegin, updateEAreaEnd;
        // Index space being updated in non-PML area.
        Int3 internalBAreaBegin, internalBAreaEnd;
        Int3 internalEAreaBegin, internalEAreaEnd;

        FP globalTime;
        FP dt;
    };

    template<GridTypes gridType>
    inline void FieldSolver<gridType>::save(std::ostream& ostr)
    {
        ostr.write((char*)&globalTime, sizeof(globalTime));
        // if (pml.get()) pml->save(ostr);  // virtual
        // if (generator.get()) generator->save(ostr);  // virtual
    }

    template<GridTypes gridType>
    inline void FieldSolver<gridType>::load(std::istream& istr)
    {
        // dt is not loaded, we get it with FieldSolver constructor
        // timeShiftE, timeShiftB, timeShiftJ are got with FieldSolver constructor too
        istr.read((char*)&globalTime, sizeof(globalTime));
        // all the next modules have already created in FieldSolver constructor
        // we need to load them only
        // if (pml.get()) pml->load(istr);  // virtual
        // if (generator.get()) generator->load(istr);  // virtual
    }

    template<GridTypes gridType>
    inline void FieldSolver<gridType>::updateDims()
    {
        updateEAreaBegin = Int3(0, 0, 0);
        updateEAreaEnd = grid->numCells;
        updateBAreaBegin = Int3(0, 0, 0);
        updateBAreaEnd = grid->numCells;
    }

    template<GridTypes gridType>
    inline void FieldSolver<gridType>::updateInternalDims()
    {
        if (pml)
        {
            for (int d = 0; d < 3; ++d)
            {
                internalBAreaBegin[d] = std::max(updateBAreaBegin[d], pml->leftPmlBorder[d]);
                internalBAreaEnd[d] = std::min(updateBAreaEnd[d], pml->rightPmlBorder[d]);
                internalEAreaBegin[d] = std::max(updateEAreaBegin[d], pml->leftPmlBorder[d]);
                internalEAreaEnd[d] = std::min(updateEAreaEnd[d], pml->rightPmlBorder[d]);
            }
        }
        else
        {
            internalBAreaBegin = updateBAreaBegin;
            internalBAreaEnd = updateBAreaEnd;
            internalEAreaBegin = updateEAreaBegin;
            internalEAreaEnd = updateEAreaEnd;
        }
    }


    template<GridTypes gridType>
    class RealFieldSolver : public FieldSolver<gridType>
    {
    public:
        RealFieldSolver(Grid<FP, gridType>* _grid, FP dt) :
            FieldSolver<gridType>(_grid, dt)
        {}

    private:
        // Copy and assignment are disallowed.
        RealFieldSolver(const RealFieldSolver&);
        RealFieldSolver& operator =(const RealFieldSolver&);
    };


    template<GridTypes gridType>
    class SpectralFieldSolver : public FieldSolver<gridType>
    {
    public:
        SpectralFieldSolver(Grid<FP, gridType>* _grid, FP dt) :
            FieldSolver<gridType>(_grid, dt)
        {
            complexGrid.reset(new SpectralGrid<FP, complexFP>(
                fourier_transform::getSizeOfComplexArray(_grid->numCells),
                _grid->globalGridDims, _grid));
            fourierTransform.initialize<gridType>(_grid, complexGrid.get());
        }

        void doFourierTransformB(fourier_transform::Direction direction);
        void doFourierTransformE(fourier_transform::Direction direction);
        void doFourierTransformJ(fourier_transform::Direction direction);
        void doFourierTransform(fourier_transform::Direction direction);

        FP3 getWaveVector(const Int3& ind);

        void updateDims();

        std::unique_ptr<SpectralGrid<FP, complexFP>> complexGrid;

        Int3 updateComplexEAreaBegin, updateComplexEAreaEnd;
        Int3 updateComplexBAreaBegin, updateComplexBAreaEnd;

        FourierTransformGrid fourierTransform;

    private:
        // Copy and assignment are disallowed.
        SpectralFieldSolver(const SpectralFieldSolver&);
        SpectralFieldSolver& operator=(const SpectralFieldSolver&);
    };

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransformB(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::x, direction);
        fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::y, direction);
        fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::z, direction);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransformE(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(FieldEnum::E, CoordinateEnum::x, direction);
        fourierTransform.doFourierTransform(FieldEnum::E, CoordinateEnum::y, direction);
        fourierTransform.doFourierTransform(FieldEnum::E, CoordinateEnum::z, direction);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransformJ(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(FieldEnum::J, CoordinateEnum::x, direction);
        fourierTransform.doFourierTransform(FieldEnum::J, CoordinateEnum::y, direction);
        fourierTransform.doFourierTransform(FieldEnum::J, CoordinateEnum::z, direction);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransform(fourier_transform::Direction direction)
    {
        doFourierTransformE(direction);
        doFourierTransformB(direction);
        doFourierTransformJ(direction);
    }

    template<GridTypes gridType>
    inline FP3 SpectralFieldSolver<gridType>::getWaveVector(const Int3& ind)
    {
        FP kx = (2 * constants::pi * ((ind.x <= this->grid->numCells.x / 2) ? ind.x : ind.x - this->grid->numCells.x)) /
            (this->grid->steps.x * this->grid->numCells.x);
        FP ky = (2 * constants::pi * ((ind.y <= this->grid->numCells.y / 2) ? ind.y : ind.y - this->grid->numCells.y)) /
            (this->grid->steps.y * this->grid->numCells.y);
        FP kz = (2 * constants::pi * ((ind.z <= this->grid->numCells.z / 2) ? ind.z : ind.z - this->grid->numCells.z)) /
            (this->grid->steps.z * this->grid->numCells.z);
        return FP3(kx, ky, kz);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::updateDims()
    {
        this->updateEAreaBegin = Int3(0, 0, 0);
        this->updateEAreaEnd = this->grid->numCells;
        this->updateBAreaBegin = Int3(0, 0, 0);
        this->updateBAreaEnd = this->grid->numCells;
        this->updateComplexEAreaBegin = Int3(0, 0, 0);
        this->updateComplexEAreaEnd = this->complexGrid->numCells;
        this->updateComplexBAreaBegin = Int3(0, 0, 0);
        this->updateComplexBAreaEnd = this->complexGrid->numCells;
    }

};