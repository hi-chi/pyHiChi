#pragma once
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"

namespace pfc
{
	template<GridTypes gridTypes>
	class FieldSolver;
	template<GridTypes gridTypes>
	class RealFieldSolver;

	template<GridTypes gridTypes>
	class FieldGenerator
	{
	public:

		FieldGenerator(RealFieldSolver<gridTypes>* fieldSolver);

        // copy constructor, other fieldSolver is possible
        FieldGenerator(const FieldGenerator& gen, RealFieldSolver<gridTypes>* fieldSolver = 0);

		virtual void generateB();
		virtual void generateE();
		RealFieldSolver<gridTypes>* fieldSolver;

        virtual FieldGenerator* createInstance(RealFieldSolver<gridTypes>* fieldSolver) = 0;

	private:

		// major index is index of edge, minor index is index of component
		std::function<FP(FP, FP, FP, FP)> eLeft[3][3];
		std::function<FP(FP, FP, FP, FP)> eRight[3][3];
		std::function<FP(FP, FP, FP, FP)> bLeft[3][3];
		std::function<FP(FP, FP, FP, FP)> bRight[3][3];
		FP3 leftCoeff;
		FP3 rightCoeff;
	};

	typedef FieldGenerator<GridTypes::YeeGridType> FieldGeneratorYee;

	template<GridTypes gridTypes>
	FieldGenerator<gridTypes>::FieldGenerator(RealFieldSolver<gridTypes>* _fieldSolver)
	{
		fieldSolver = _fieldSolver;
		leftCoeff = FP3(0, 0, 0);
		rightCoeff = FP3(0, 0, 0);
		for (int d = 0; d < 3; ++d)
		{
			//in seq if fieldGeneration in area
			leftCoeff[d] = 1;
			rightCoeff[d] = 1;
		}

		fieldSolver->generator.reset(this);
	}

    template<GridTypes gridTypes>
    inline FieldGenerator<gridTypes>::FieldGenerator(const FieldGenerator & gen,
        RealFieldSolver<gridTypes>* fieldSolver = 0)
    {
        if (fieldSolver)
            this->fieldSolver = fieldSolver;
        else this->fieldSolver = gen.fieldSolver;

        leftCoeff = gen.leftCoeff;
        rightCoeff = gen.rightCoeff;

        for (int f = 0; f < 3; ++f)
            for (int d = 0; d < 3; ++d) {
                eLeft[f][d] = gen.eLeft[f][d];
                eRight[f][d] = gen.eRight[f][d];
                bLeft[f][d] = gen.bLeft[f][d];
                bRight[f][d] = gen.bRight[f][d];
            }

        this->fieldSolver->generator.reset(this);
    }

	template<GridTypes gridTypes>
	void FieldGenerator<gridTypes>::generateB()
	{
		Grid<FP, gridTypes>* grid = fieldSolver->grid;
		const FP time = fieldSolver->time + grid->timeShiftB;
		const FP cdt = constants::c * grid->dt;
		const FP3 norm_coeffs = FP3(cdt, cdt, cdt) / grid->steps;
		for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
		{
			int dim1 = (dim0 + 1) % 3;
			int dim2 = (dim0 + 2) % 3;
			int begin1 = fieldSolver->internalBAreaBegin[dim1];
			int begin2 = fieldSolver->internalBAreaBegin[dim2];
			int end1 = fieldSolver->internalBAreaEnd[dim1];
			int end2 = fieldSolver->internalBAreaEnd[dim2];
//#pragma omp parallel for collapse(2)
			for (int j = begin1; j < end1; j++)
				for (int k = begin2; k < end2; k++)
				{
					// Adjust indexes for symmetry of generation coordinates
					Int3 index;
					index[dim0] = fieldSolver->internalBAreaBegin[dim0] + 1;
					index[dim1] = j;
					index[dim2] = k;
					Int3 indexes[3] = { index, index, index };
					indexes[dim1][dim0]++;
					indexes[dim2][dim0]++;
					FP3 bxCoords = grid->BxPosition(indexes[0].x, indexes[0].y,
						indexes[0].z);
					FP3 byCoords = grid->ByPosition(indexes[1].x, indexes[1].y,
						indexes[1].z);
					FP3 bzCoords = grid->BzPosition(indexes[2].x, indexes[2].y,
						indexes[2].z);
					FP coeff = leftCoeff[dim0] * norm_coeffs[dim0];
					grid->Bx(indexes[0]) += coeff * bLeft[dim0][0](
						bxCoords.x, bxCoords.y, bxCoords.z, time);
					grid->By(indexes[1]) += coeff * bLeft[dim0][1](
						byCoords.x, byCoords.y, byCoords.z, time);
					grid->Bz(indexes[2]) += coeff * bLeft[dim0][2](
						bzCoords.x, bzCoords.y, bzCoords.z, time);

					index[dim0] = fieldSolver->internalBAreaEnd[dim0] - 2;
					bxCoords = grid->BxPosition(index.x, index.y, index.z);
					byCoords = grid->ByPosition(index.x, index.y, index.z);
					bzCoords = grid->BzPosition(index.x, index.y, index.z);
					coeff = rightCoeff[dim0] * norm_coeffs[dim0];
					grid->Bx(index) += coeff * bRight[dim0][0](
						bxCoords.x, bxCoords.y, bxCoords.z, time);
					grid->By(index) += coeff * bRight[dim0][1](
						byCoords.x, byCoords.y, byCoords.z, time);
					grid->Bz(index) += coeff * bRight[dim0][2](
						bzCoords.x, bzCoords.y, bzCoords.z, time);
				}
		}
	}


	template<GridTypes gridTypes>
	void FieldGenerator<gridTypes>::generateE()
	{
		Grid<FP, gridTypes> * grid = fieldSolver->grid;
		const FP time = fieldSolver->time + grid->timeShiftB;
		const FP cdt = constants::c * grid->dt;
		const FP3 norm_coeffs = FP3(cdt, cdt, cdt) / grid->steps;
		for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
		{
			int dim1 = (dim0 + 1) % 3;
			int dim2 = (dim0 + 2) % 3;
			int begin1 = fieldSolver->internalEAreaBegin[dim1];
			int begin2 = fieldSolver->internalEAreaBegin[dim2];
			int end1 = fieldSolver->internalEAreaEnd[dim1];
			int end2 = fieldSolver->internalEAreaEnd[dim2];
//#pragma omp parallel for collapse(2)
			for (int j = begin1; j < end1; j++)
				for (int k = begin2; k < end2; k++)
				{
					// Adjust indexes for symmetry of generation coordinates, use B indexes for consistency
					/// Видимо, не нужно добавлять здесь 1, а ниже вычитать 1 вместо 2
					Int3 index;
					index[dim0] = fieldSolver->internalBAreaBegin[dim0] + 1;
					index[dim1] = j;
					index[dim2] = k;
					Int3 indexes[3] = { index, index, index };
					indexes[dim0][dim0]++;
					FP3 exCoords = grid->ExPosition(indexes[0].x, indexes[0].y,
						indexes[0].z);
					FP3 eyCoords = grid->EyPosition(indexes[1].x, indexes[1].y,
						indexes[1].z);
					FP3 ezCoords = grid->EzPosition(indexes[2].x, indexes[2].y,
						indexes[2].z);
					FP coeff = leftCoeff[dim0] * norm_coeffs[dim0];
					grid->Ex(indexes[0]) += coeff * eLeft[dim0][0](
						exCoords.x, exCoords.y, exCoords.z, time);
					grid->Ey(indexes[1]) += coeff * eLeft[dim0][1](
						eyCoords.x, eyCoords.y, eyCoords.z, time);
					grid->Ez(indexes[2]) += coeff * eLeft[dim0][2](
						ezCoords.x, ezCoords.y, ezCoords.z, time);

					index[dim0] = fieldSolver->internalBAreaEnd[dim0] - 2;
					exCoords = grid->ExPosition(index.x, index.y, index.z);
					eyCoords = grid->EyPosition(index.x, index.y, index.z);
					ezCoords = grid->EzPosition(index.x, index.y, index.z);
					coeff = rightCoeff[dim0] * norm_coeffs[dim0];
					grid->Ex(index) += coeff * eRight[dim0][0](
						exCoords.x, exCoords.y, exCoords.z, time);
					grid->Ey(index) += coeff * eRight[dim0][1](
						eyCoords.x, eyCoords.y, eyCoords.z, time);
					grid->Ez(index) += coeff * eRight[dim0][2](
						ezCoords.x, ezCoords.y, ezCoords.z, time);
				}
		}
	}


	template<GridTypes gridTypes>
	class PeriodicalFieldGenerator : public FieldGenerator<gridTypes>
	{
	public:
		PeriodicalFieldGenerator(RealFieldSolver<gridTypes>* fieldSolver) :
			FieldGenerator<gridTypes>(fieldSolver) {
		}

        // copy constructor, other fieldSolver is possible
        PeriodicalFieldGenerator(const PeriodicalFieldGenerator& gen, RealFieldSolver<gridTypes>* fieldSolver = 0) :
            FieldGenerator<gridTypes>(gen, fieldSolver) {
        }

		virtual void generateB();
		virtual void generateE();

        FieldGenerator* createInstance(RealFieldSolver<gridTypes>* fieldSolver) override {
            return new PeriodicalFieldGenerator(*this, fieldSolver);
        }
	};

	typedef PeriodicalFieldGenerator<GridTypes::YeeGridType> PeriodicalFieldGeneratorYee;

	template<GridTypes gridTypes>
	void PeriodicalFieldGenerator<gridTypes>::generateB()
	{
		Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
		for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
		{
			int dim1 = (dim0 + 1) % 3;
			int dim2 = (dim0 + 2) % 3;
			int begin1 = this->fieldSolver->internalEAreaBegin[dim1];
			int begin2 = this->fieldSolver->internalEAreaBegin[dim2];
			int end1 = this->fieldSolver->internalEAreaEnd[dim1];
			int end2 = this->fieldSolver->internalEAreaEnd[dim2];
//#pragma omp parallel for collapse(2)
			for (int j = begin1; j < end1; j++)
				for (int k = begin2; k < end2; k++)
				{
					// Adjust indexes for symmetry of generation coordinates
					Int3 indexL, indexR;
					indexL[dim0] = this->fieldSolver->internalEAreaBegin[dim0];
					indexL[dim1] = j;
					indexL[dim2] = k;
					indexR[dim0] = this->fieldSolver->internalEAreaEnd[dim0] - 3;
					indexR[dim1] = j;
					indexR[dim2] = k;


					grid->Bx(indexL) = grid->Bx(indexR);
					grid->By(indexL) = grid->By(indexR);
					grid->Bz(indexL) = grid->Bz(indexR);

					indexL[dim0]++;
					indexR[dim0]++;

					grid->Bx(indexR) = grid->Bx(indexL);
					grid->By(indexR) = grid->By(indexL);
					grid->Bz(indexR) = grid->Bz(indexL);
				}
		}
	}


	template<GridTypes gridTypes>
	void PeriodicalFieldGenerator<gridTypes>::generateE()
	{
		Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
		for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
		{
			int dim1 = (dim0 + 1) % 3;
			int dim2 = (dim0 + 2) % 3;
			int begin1 = 0;
			int begin2 = 0;
			int end1 = grid->numCells[dim1];
			int end2 = grid->numCells[dim2];
//#pragma omp parallel for collapse(2)
			for (int j = begin1; j < end1; j++)
				for (int k = begin2; k < end2; k++)
				{
					// Adjust indexes for symmetry of generation coordinates
					Int3 indexL, indexR;
					indexL[dim0] = this->fieldSolver->internalEAreaBegin[dim0] + 1;
					indexL[dim1] = j;
					indexL[dim2] = k;
					indexR[dim0] = this->fieldSolver->internalEAreaEnd[dim0] - 2;
					indexR[dim1] = j;
					indexR[dim2] = k;

					grid->Ex(indexL) = grid->Ex(indexR);
					grid->Ey(indexL) = grid->Ey(indexR);
					grid->Ez(indexL) = grid->Ez(indexR);

					indexL[dim0]++;
					indexR[dim0]++;

					grid->Ex(indexR) = grid->Ex(indexL);
					grid->Ey(indexR) = grid->Ey(indexL);
					grid->Ez(indexR) = grid->Ez(indexL);
				}
		}
	}

	template<GridTypes gridTypes>
    class ReflectFieldGenerator : public FieldGenerator<gridTypes>
    {
    public:
        ReflectFieldGenerator(RealFieldSolver<gridTypes>* fieldSolver) :
            FieldGenerator<gridTypes>(fieldSolver) {
        };

        // copy constructor, other fieldSolver is possible
        ReflectFieldGenerator(const ReflectFieldGenerator& gen, RealFieldSolver<gridTypes>* fieldSolver = 0) :
            FieldGenerator<gridTypes>(gen, fieldSolver) {
        }

        virtual void generateB();
        virtual void generateE();

        FieldGenerator* createInstance(RealFieldSolver<gridTypes>* fieldSolver) override {
            return new ReflectFieldGenerator(*this, fieldSolver);
        }
    };

	typedef ReflectFieldGenerator<GridTypes::YeeGridType> ReflectFieldGeneratorYee;

	template<GridTypes gridTypes>
	void ReflectFieldGenerator<gridTypes>::generateB()
	{
		Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
		for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
		{
			int dim1 = (dim0 + 1) % 3;
			int dim2 = (dim0 + 2) % 3;
			int begin1 = this->fieldSolver->internalBAreaBegin[dim1];
			int begin2 = this->fieldSolver->internalBAreaBegin[dim2];
			int end1 = this->fieldSolver->internalBAreaEnd[dim1];
			int end2 = this->fieldSolver->internalBAreaEnd[dim2];
//#pragma omp parallel for collapse(2)
			for (int j = begin1; j < end1; j++)
				for (int k = begin2; k < end2; k++)
				{
					// Adjust indexes for symmetry of generation coordinates
					Int3 indexL, indexR;
					indexL[dim0] = this->fieldSolver->internalBAreaBegin[dim0];
					indexL[dim1] = j;
					indexL[dim2] = k;
					indexR[dim0] = this->fieldSolver->internalBAreaEnd[dim0] - 1;
					indexR[dim1] = j;
					indexR[dim2] = k;

					grid->Bx(indexR) = (FP)0.0;
					grid->By(indexR) = (FP)0.0;
					grid->Bz(indexR) = (FP)0.0;
					grid->Bx(indexL) = (FP)0.0;
					grid->By(indexL) = (FP)0.0;
					grid->Bz(indexL) = (FP)0.0;
				}
		}
	}


	template<GridTypes gridTypes>
	void ReflectFieldGenerator<gridTypes>::generateE()
	{
		Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
		for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
		{
			int dim1 = (dim0 + 1) % 3;
			int dim2 = (dim0 + 2) % 3;
			int begin1 = this->fieldSolver->internalEAreaBegin[dim1];
			int begin2 = this->fieldSolver->internalEAreaBegin[dim2];
			int end1 = this->fieldSolver->internalEAreaEnd[dim1];
			int end2 = this->fieldSolver->internalEAreaEnd[dim2];
//#pragma omp parallel for collapse(2)
			for (int j = begin1; j < end1; j++)
				for (int k = begin2; k < end2; k++)
				{
					// Adjust indexes for symmetry of generation coordinates
					Int3 indexL, indexR;
					indexL[dim0] = this->fieldSolver->internalEAreaBegin[dim0];
					indexL[dim1] = j;
					indexL[dim2] = k;
					indexR[dim0] = this->fieldSolver->internalEAreaEnd[dim0] - 1;
					indexR[dim1] = j;
					indexR[dim2] = k;


					grid->Ex(indexL) = (FP)0.0;
					grid->Ey(indexL) = (FP)0.0;
					grid->Ez(indexL) = (FP)0.0;
					grid->Ex(indexR) = (FP)0.0;
					grid->Ey(indexR) = (FP)0.0;
					grid->Ez(indexR) = (FP)0.0;
				}
		}
	}
}