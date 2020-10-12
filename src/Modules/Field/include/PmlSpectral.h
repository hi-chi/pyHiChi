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
            return static_cast<SpectralFieldSolver<gridTypes>*>(this->fieldSolver);
        }

        void computeSplitField(ScalarField<complexFP>& field,
            std::vector<FP>& splitField, std::vector<FP>& coeffs, Coordinate splitCoord);

        FP getZeroComponentCoeff(Coordinate propDir);
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
        coeffx.resize(this->numCells);
        coeffy.resize(this->numCells);
        coeffz.resize(this->numCells);
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
                    if (i == 0 && j == 0 && k == 0) {
                        this->tmpFieldComplex(0, 0, 0) = this->getZeroComponentCoeff(splitCoord) * field(0, 0, 0);
                    }
                    else {
                        FP3 K = SpectralFieldSolver<gridTypes>::getWaveVector(Int3(i, j, k), gridSize, gridStep);
                        FP coeff = K[splitCoord] * K[splitCoord] / K.norm2();
                        this->tmpFieldComplex(i, j, k) = coeff * field(i, j, k);
                    }
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
        Grid<FP, gridTypes> * grid = this->getFieldSolver()->grid;
        for (int idx = 0; idx < this->numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            grid->Bx(i, j, k) = this->bxy[idx] + this->bxz[idx];
            grid->By(i, j, k) = this->byz[idx] + this->byx[idx];
            grid->Bz(i, j, k) = this->bzx[idx] + this->bzy[idx];

            grid->Ex(i, j, k) = this->exy[idx] + this->exz[idx];
            grid->Ey(i, j, k) = this->eyz[idx] + this->eyx[idx];
            grid->Ez(i, j, k) = this->ezx[idx] + this->ezy[idx];
        }
    }

    template<GridTypes gridTypes>
    inline FP PmlSpectral<gridTypes>::getZeroComponentCoeff(Coordinate propDir)
    {
        if (propDir == Coordinate::x) return 1.0;
        //Grid<complexFP, gridTypes>* grid = this->getFieldSolver()->complexGrid;
        //const FP eps = 1e-30;  // only for double
        //// zero component is always real
        //FP Ex = grid->Ex(0, 0, 0).real, Ey = grid->Ey(0, 0, 0).real, Ez = grid->Ez(0, 0, 0).real;
        //FP Bx = grid->Bx(0, 0, 0).real, By = grid->By(0, 0, 0).real, Bz = grid->Bz(0, 0, 0).real;
        //{
        //    FP ExBy_BxEy = Ex * By - Bx * Ey;
        //    if (ExBy_BxEy > eps || ExBy_BxEy < -eps) {
        //        FP kx = (Bx * Bz + Ex * Ez) / ExBy_BxEy;
        //        FP ky = (By * Bz + Ey * Ez) / ExBy_BxEy;
        //        switch (propDir) {
        //        case Coordinate::x:
        //            return kx * kx;
        //        case Coordinate::y:
        //            return ky * ky;
        //        default:  // Coordinate::z
        //            return 1.0 - kx * kx - ky * ky;
        //        }
        //    }
        //}
        //{
        //    FP EzBx_BzEx = Ez * Bx - Bz * Ex;
        //    if (EzBx_BzEx > eps || EzBx_BzEx < -eps) {
        //        FP kx = (Bx * By + Ex * Ey) / EzBx_BzEx;
        //        FP kz = (By * Bz + Ey * Ez) / EzBx_BzEx;
        //        switch (propDir) {
        //        case Coordinate::x:
        //            return kx * kx;
        //        case Coordinate::z:
        //            return kz * kz;
        //        default:  // Coordinate::y
        //            return 1.0 - kx * kx - kz * kz;
        //        }
        //    }
        //}
        //{
        //    FP EyBz_ByEz = Ey * Bz - By * Ez;
        //    if (EyBz_ByEz > eps || EyBz_ByEz < -eps) {
        //        FP ky = (Bx * By + Ex * Ey) / EyBz_ByEz;
        //        FP kz = (Bx * Bz + Ex * Ez) / EyBz_ByEz;
        //        switch (propDir) {
        //        case Coordinate::y:
        //            return ky * ky;
        //        case Coordinate::z:
        //            return kz * kz;
        //        default:  // Coordinate::x
        //            return 1.0 - ky * ky - kz * kz;
        //        }
        //    }
        //}
        return 0.0;
    }
}