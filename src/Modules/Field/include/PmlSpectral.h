#pragma once
#include "Grid.h"
#include "FieldSolver.h"
#include "Pml.h"
#include "Constants.h"

namespace pfc {

    template<GridTypes gridTypes>
    class PmlSpectral : public Pml<gridTypes>
    {
    public:

        PmlSpectral(SpectralFieldSolver<gridTypes>* fieldSolver, Int3 sizePML);

        void computeCoeffs();

        void updateSplit();

        void sumSplit();

        int numCells; // total number of PML cells
        std::vector<Int3> cellIndex; // natural 3d indexes of cells in PML
        std::vector<FP> coeffx, coeffy, coeffz;

        ScalarField<FP> tmpFieldReal;
        ScalarField<complexFP> tmpFieldComplex;  // has shared memory with tmpFieldReal
        FourierTransformField fourierTransform;

    protected:
        SpectralFieldSolver<gridTypes>* getFieldSolver() {
            return (SpectralFieldSolver<gridTypes>*)this->fieldSolver;
        }

        void computeSplitField(ScalarField<complexFP>& field,
            std::vector<FP>& splitField, std::vector<FP>& coeffs, Coordinate splitCoord);
    };

    template<GridTypes gridTypes>
    inline PmlSpectral<gridTypes>::PmlSpectral(SpectralFieldSolver<gridTypes>* _fieldSolver, Int3 _sizePML) :
        Pml<gridTypes>((FieldSolver<gridTypes>*) _fieldSolver, _sizePML),
        tmpFieldReal(fieldSolver->grid->sizeStorage),
        tmpFieldComplex(reinterpret_cast<complexFP*>(tmpFieldReal.getData()),
            static_cast<SpectralFieldSolver<gridTypes>*>(fieldSolver)->complexGrid->sizeStorage)
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;

        const Int3 leftPmlEnd = this->leftDims;
        const Int3 rightPmlBegin = grid->numCells - this->rightDims;
        const int maxSize = grid->numCells.volume() -
            (grid->numCells - this->leftDims - this->rightDims).volume();

        Int3 begin = this->fieldSolver->updateEAreaBegin; // num nodes B = num nodes E
        Int3 end = this->fieldSolver->updateEAreaEnd;
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
                    bool xBoundaryPml = (i < leftPmlEnd.x) || (i >= rightPmlBegin.x);
                    bool yBoundaryPml = (j < leftPmlEnd.y) || (j >= rightPmlBegin.y);
                    bool zBoundaryPml = (k < leftPmlEnd.z) || (k >= rightPmlBegin.z);
                    if (xBoundaryPml || yBoundaryPml || zBoundaryPml)
                    {
                        cellIndex.push_back(Int3(i, j, k));
                        int index = k + (j + i * grid->numCells.y) * grid->numCells.z;
                    }
                }

        numCells = (int)cellIndex.size();
        this->initializeSplitFieldsE(numCells);
        this->initializeSplitFieldsB(numCells);

        computeCoeffs();

        fourierTransform.initialize(&tmpFieldReal, &tmpFieldComplex, this->fieldSolver->grid->numCells);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::computeCoeffs()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        coeffx.resize(numCells);
        coeffy.resize(numCells);
        coeffz.resize(numCells);
        for (int idx = 0; idx < numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            FP3 coords = grid->ExPosition(i, j, k);
            FP3 sigma = this->computeSigma(coords);

            coeffx[idx] = exp(-sigma.x * this->fieldSolver->dt);
            coeffy[idx] = exp(-sigma.y * this->fieldSolver->dt);
            coeffz[idx] = exp(-sigma.z * this->fieldSolver->dt);
        }
    }

    // TODO: vectorisation
    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::computeSplitField(ScalarField<complexFP>& field,
        std::vector<FP>& splitField, std::vector<FP>& coeffs, Coordinate splitCoord)
    {
        SpectralFieldSolver<gridTypes>* fs = this->getFieldSolver();
        Int3 gridSize = fs->grid->numCells;
        Int3 complexGridSize = fs->complexGrid->numCells;
        FP3 gridStep = fs->grid->steps;

        for (int i = 0; i < complexGridSize.x; i++)
            for (int j = 0; j < complexGridSize.y; j++)
                for (int k = 0; k < complexGridSize.z; k++) {
                    FP3 K = SpectralFieldSolver<gridTypes>::getWaveVector(Int3(i, j, k), gridSize, gridStep);
                    FP coeff = K[splitCoord] * K[splitCoord] / K.norm2();
                    this->tmpFieldComplex(i, j, k) = coeff * field(i, j, k);
                }

        this->fourierTransform.doInverseFourierTransform();

        for (int idx = 0; idx < this->numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            splitField[idx] = coeffs[idx] * this->tmpFieldReal(i, j, k);
        }
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::updateSplit()
    {
        this->computeSplitField(this->getFieldSolver()->complexGrid->Bx, this->bxy,
            this->coeffy, Coordinate::y);
        this->computeSplitField(this->getFieldSolver()->complexGrid->Bx, this->bxz,
            this->coeffz, Coordinate::z);

        this->computeSplitField(this->getFieldSolver()->complexGrid->By, this->byz,
            this->coeffz, Coordinate::z);
        this->computeSplitField(this->getFieldSolver()->complexGrid->By, this->byx,
            this->coeffx, Coordinate::x);

        this->computeSplitField(this->getFieldSolver()->complexGrid->Bz, this->bzx,
            this->coeffx, Coordinate::x);
        this->computeSplitField(this->getFieldSolver()->complexGrid->Bz, this->bzy,
            this->coeffy, Coordinate::y);

        this->computeSplitField(this->getFieldSolver()->complexGrid->Ex, this->exy,
            this->coeffy, Coordinate::y);
        this->computeSplitField(this->getFieldSolver()->complexGrid->Ex, this->exz,
            this->coeffz, Coordinate::z);

        this->computeSplitField(this->getFieldSolver()->complexGrid->Ey, this->eyz,
            this->coeffz, Coordinate::z);
        this->computeSplitField(this->getFieldSolver()->complexGrid->Ey, this->eyx,
            this->coeffx, Coordinate::x);

        this->computeSplitField(this->getFieldSolver()->complexGrid->Ez, this->ezx,
            this->coeffx, Coordinate::x);
        this->computeSplitField(this->getFieldSolver()->complexGrid->Ez, this->ezy,
            this->coeffy, Coordinate::y);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumSplit()
    {
        for (int idx = 0; idx < this->numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            this->getFieldSolver()->grid->Bx(i, j, k) = this->bxy[idx] + this->bxz[idx];
            this->getFieldSolver()->grid->By(i, j, k) = this->byz[idx] + this->byx[idx];
            this->getFieldSolver()->grid->Bz(i, j, k) = this->bzx[idx] + this->bzy[idx];

            this->getFieldSolver()->grid->Ex(i, j, k) = this->exy[idx] + this->exz[idx];
            this->getFieldSolver()->grid->Ey(i, j, k) = this->eyz[idx] + this->eyx[idx];
            this->getFieldSolver()->grid->Ez(i, j, k) = this->ezx[idx] + this->ezy[idx];
        }
    }
}