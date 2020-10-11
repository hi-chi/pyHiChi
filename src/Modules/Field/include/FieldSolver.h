#pragma once
#include "Grid.h"
#include "FieldGenerator.h"
#include "Vectors.h"
#include "FourierTransform.h"

#include <algorithm>
#include <memory>

namespace pfc {
    template<GridTypes gridType>
    class Pml;

    // Base class for field solvers on Grid.
    // The main method doStep uses template method pattern.
    template<GridTypes gridType>
    class FieldSolver
    {

    public:

        FieldSolver(Grid<FP, gridType>* _grid, FP dt,
            FP timeShiftE, FP timeShiftB, FP timeShiftJ):
            dt(dt),
            timeShiftE(timeShiftE),
            timeShiftB(timeShiftB),
            timeShiftJ(timeShiftJ)
        {
            this->grid = _grid;
            this->globalTime = (FP)0.0;
            this->generator.reset(nullptr);
        }

        /*void setPML(int sizePMLx, int sizePMLy, int sizePMLz)
        {
            pml.reset(new Pml(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
        }*/

        Grid<FP, gridType> * grid;

        std::unique_ptr<Pml<gridType>> pml;
        std::unique_ptr<FieldGenerator<gridType>> generator;

        // Index space being updated in form [begin, end).
        Int3 updateBAreaBegin, updateBAreaEnd;
        Int3 updateEAreaBegin, updateEAreaEnd;
        // Index space being updated in non-PML area.
        Int3 internalBAreaBegin, internalBAreaEnd;
        Int3 internalEAreaBegin, internalEAreaEnd;

        void updateDims();
        void updateInternalDims();

        FP globalTime;
        FP dt;
        // difference between E and B
        FP timeShiftE, timeShiftB, timeShiftJ;
    };
    

    template<GridTypes gridType>
    inline void FieldSolver<gridType>::updateDims()
    {
        updateEAreaBegin = Int3(0, 0, 0);
        updateEAreaEnd = grid->numCells - grid->getNumExternalRightCells() / 2;
        updateBAreaBegin = grid->getNumExternalLeftCells();
        updateBAreaEnd = grid->numCells - updateBAreaBegin;
    }

    template<GridTypes gridType>
    inline void FieldSolver<gridType>::updateInternalDims()
    {
        if (pml.get())
        {
            for (int d = 0; d < 3; ++d)
            {
                internalBAreaBegin[d] = std::max(updateBAreaBegin[d], pml->leftDims[d]);
                internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
                internalEAreaBegin[d] = std::max(updateEAreaBegin[d], pml->leftDims[d]);
                internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
            }
        }
        else
        {
            for (int d = 0; d < 3; ++d)
            {
                internalBAreaBegin[d] = updateBAreaBegin[d];
                internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
                    grid->numCells[d]);
                internalEAreaBegin[d] = updateEAreaBegin[d];
                internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
                    grid->numCells[d]);
            }
        }
    }


    template<GridTypes gridType>
    class RealFieldSolver : public FieldSolver<gridType>
    {
    public:
        RealFieldSolver(Grid<FP, gridType>* _grid, FP dt,
            FP timeShiftE, FP timeShiftB, FP timeShiftJ):
            FieldSolver<gridType>(_grid, dt, timeShiftE, timeShiftB, timeShiftJ)
        {}

    private:
        // Copy and assignment are disallowed.
        RealFieldSolver(const RealFieldSolver &);
        RealFieldSolver & operator =(const RealFieldSolver &);
    };


    template<GridTypes gridType>
    class SpectralFieldSolver : public FieldSolver<gridType>
    {
    public:
        SpectralFieldSolver(Grid<FP, gridType>* _grid, FP dt,
            FP timeShiftE, FP timeShiftB, FP timeShiftJ) :
            FieldSolver<gridType>(_grid, dt, timeShiftE, timeShiftB, timeShiftJ)
        {
            complexGrid = new Grid<complexFP, gridType>(fourier_transform::getSizeOfComplexArray(_grid->numCells),
                fourier_transform::getSizeOfComplexArray(_grid->globalGridDims), _grid);
            fourierTransform.initialize<gridType>(_grid, complexGrid);
        }

        ~SpectralFieldSolver() {
            delete complexGrid;
        }

        void doFourierTransformB(fourier_transform::Direction direction);
        void doFourierTransformE(fourier_transform::Direction direction);
        void doFourierTransformJ(fourier_transform::Direction direction);
        void doFourierTransform(fourier_transform::Direction direction);

        static FP getWaveVectorComponent(int ind, int numCells, FP gridStep);
        static FP3 getWaveVector(const Int3 & ind, const Int3 & numCells, const FP3 & gridStep);
        FP3 getWaveVector(const Int3 & ind);

        void updateDims();

        Grid<complexFP, gridType> * complexGrid;

        Int3 updateComplexEAreaBegin, updateComplexEAreaEnd;
        Int3 updateComplexBAreaBegin, updateComplexBAreaEnd;

        FourierTransformGrid fourierTransform;

    private:
        // Copy and assignment are disallowed.
        SpectralFieldSolver(const SpectralFieldSolver &);
        SpectralFieldSolver & operator=(const SpectralFieldSolver &);
    };

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransformB(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(B, x, direction);
        fourierTransform.doFourierTransform(B, y, direction);
        fourierTransform.doFourierTransform(B, z, direction);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransformE(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(E, x, direction);
        fourierTransform.doFourierTransform(E, y, direction);
        fourierTransform.doFourierTransform(E, z, direction);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransformJ(fourier_transform::Direction direction)
    {
        fourierTransform.doFourierTransform(J, x, direction);
        fourierTransform.doFourierTransform(J, y, direction);
        fourierTransform.doFourierTransform(J, z, direction);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransform(fourier_transform::Direction direction)
    {
        doFourierTransformE(direction);
        doFourierTransformB(direction);
        doFourierTransformJ(direction);
    }

    template<GridTypes gridType>
    forceinline FP SpectralFieldSolver<gridType>::getWaveVectorComponent(int ind, int numCells, FP gridStep)
    {
        return (2 * constants::pi*((ind <= numCells / 2) ? ind : ind - numCells)) / (gridStep * numCells);
    }

    template<GridTypes gridType>
    forceinline FP3 SpectralFieldSolver<gridType>::getWaveVector(const Int3 & ind,
        const Int3 & numCells, const FP3 & steps)
    {
        return FP3(SpectralFieldSolver<gridType>::getWaveVectorComponent(ind.x, numCells.x, steps.x),
            SpectralFieldSolver<gridType>::getWaveVectorComponent(ind.y, numCells.y, steps.y),
            SpectralFieldSolver<gridType>::getWaveVectorComponent(ind.z, numCells.z, steps.z)
        );
    }

    template<GridTypes gridType>
    forceinline FP3 SpectralFieldSolver<gridType>::getWaveVector(const Int3 & ind)
    {
        return SpectralFieldSolver<gridType>::getWaveVector(ind, this->grid->numCells, this->grid->steps);
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