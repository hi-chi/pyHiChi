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

		FieldSolver(Grid<FP, gridType>* _grid)
		{
			this->grid = _grid;
			time = (FP)0.0;
			generator.reset(nullptr);
		}

		/*void setPML(int sizePMLx, int sizePMLy, int sizePMLz)
		{
			pml.reset(new Pml(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
		}*/

		Grid<FP, gridType> * grid;

		std::auto_ptr<Pml<gridType>> pml;
		std::auto_ptr<FieldGenerator<gridType>> generator;

		// Index space being updated in form [begin, end).
		Int3 updateBAreaBegin, updateBAreaEnd;
		Int3 updateEAreaBegin, updateEAreaEnd;
		// Index space being updated in non-PML area.
		Int3 internalBAreaBegin, internalBAreaEnd;
		Int3 internalEAreaBegin, internalEAreaEnd;

		void updateDims();
		void updateInternalDims();

		FP time;
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
		RealFieldSolver(Grid<FP, gridType>* _grid):
			FieldSolver<gridType>(_grid)
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
        SpectralFieldSolver(Grid<FP, gridType>* _grid) :
            FieldSolver<gridType>(_grid)
        {
            complexGrid = new Grid<complexFP, gridType>(FourierTransform::getSizeOfComplex(_grid->numCells),
                _grid->dt, FourierTransform::getSizeOfComplex(_grid->globalGridDims));
            fourierTransform.initialize<gridType>(_grid, complexGrid);
        }

        ~SpectralFieldSolver() {
            delete complexGrid;
        }

        void doFourierTransformB(FourierTransformDirection direction);
        void doFourierTransformE(FourierTransformDirection direction);
        void doFourierTransformJ(FourierTransformDirection direction);
        void doFourierTransform(FourierTransformDirection direction);

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
    inline void SpectralFieldSolver<gridType>::doFourierTransformB(FourierTransformDirection direction)
    {
        fourierTransform.doFourierTransform(B, x, direction);
        fourierTransform.doFourierTransform(B, y, direction);
        fourierTransform.doFourierTransform(B, z, direction);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransformE(FourierTransformDirection direction)
    {
        fourierTransform.doFourierTransform(E, x, direction);
        fourierTransform.doFourierTransform(E, y, direction);
        fourierTransform.doFourierTransform(E, z, direction);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransformJ(FourierTransformDirection direction)
    {
        fourierTransform.doFourierTransform(J, x, direction);
        fourierTransform.doFourierTransform(J, y, direction);
        fourierTransform.doFourierTransform(J, z, direction);
    }

    template<GridTypes gridType>
    inline void SpectralFieldSolver<gridType>::doFourierTransform(FourierTransformDirection direction)
    {
        doFourierTransformE(direction);
        doFourierTransformB(direction);
        doFourierTransformJ(direction);
    }

    template<GridTypes gridType>
    inline FP3 SpectralFieldSolver<gridType>::getWaveVector(const Int3 & ind)
    {
        FP kx = (2 * constants::pi*((ind.x <= this->grid->numCells.x / 2) ? ind.x : ind.x - this->grid->numCells.x)) /
            (this->grid->steps.x * this->grid->numCells.x);
        FP ky = (2 * constants::pi*((ind.y <= this->grid->numCells.y / 2) ? ind.y : ind.y - this->grid->numCells.y)) /
            (this->grid->steps.y * this->grid->numCells.y);
        FP kz = (2 * constants::pi*((ind.z <= this->grid->numCells.z / 2) ? ind.z : ind.z - this->grid->numCells.z)) /
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