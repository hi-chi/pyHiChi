#pragma once

#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"

#include <vector>

namespace pfc {
    template<GridTypes gridTypes>
    class FieldSolver;
    template<GridTypes gridTypes>
    class RealFieldSolver;
    template<GridTypes gridTypes>
    class SpectralFieldSolver;


    // Data for PML.
    //
    // leftDims and rightDims are number of cells, i.e. layers of b. There is one
    // layer of e on "left" outer PML boundary for which corresponding b values are
    // not used (but stored). Values of e on (all) outer PML boundary are filled
    // with zeroes, values of e in PML-internal area boundaries are copy of internal
    // area boundary e values.

    template<GridTypes gridTypes>
    class Pml
    {
    public:
        Pml(FieldSolver<gridTypes>* _fieldSolver, Int3 _sizePML);

        // only for real solvers
        virtual void updateB() {};
        virtual void updateE() {};

        Int3 sizePML;
        Int3 leftDims, rightDims;
        std::vector<FP> exy, exz, eyx, eyz, ezx, ezy; // split electric field
        std::vector<FP> bxy, bxz, byx, byz, bzx, bzy; // split magnetic field

        FP3 computeSigma(const FP3& coords);

    protected:

        FieldSolver<gridTypes> *fieldSolver;

        Int3 globalLeftDims, globalRightDims;
        FP3 maxSigma;
        Int3 leftDists, rightDists;
        FP n;

        void initializeSplitFieldsE(int sizeE);
        void initializeSplitFieldsB(int sizeB);

    };

    template<GridTypes gridTypes>
    inline Pml<gridTypes>::Pml(FieldSolver<gridTypes>* _fieldSolver, Int3 _sizePML) :
        sizePML(_sizePML), fieldSolver(_fieldSolver)
    {
        Grid<FP, gridTypes>* grid = fieldSolver->grid;
        globalLeftDims = _sizePML;

        for (int d = 0; d < 3; d++)
        {
            //printf("%x %x ", fieldSolver->grid, fieldSolver->grid->globalGridDims);
            //printf("%d, %d\n",globalLeftDims[d] * 2 , grid->globalGridDims[d]);
            if (globalLeftDims[d] * 2 > grid->globalGridDims[d])
            {
                printf("Check pml size: possibly, grid size is less than pml size");
                globalLeftDims[d] = 0;
            }
        }
        globalRightDims = globalLeftDims;
        for (int d = 0; d < 3; ++d)
        {
            int pmlLeftBorder = globalLeftDims[d];

            leftDims[d] = pmlLeftBorder;

            if (leftDims[d])
                leftDims[d] += grid->getNumExternalLeftCells()[d];

            int pmlRightBorder = grid->numInternalCells[d] - globalRightDims[d];
            rightDims[d] = globalRightDims[d] -
                (grid->globalGridDims[d] - grid->numInternalCells[d]);
            if (rightDims[d])
                rightDims[d] += grid->getNumExternalRightCells()[d];
        }

        leftDists = Int3(0, 0, 0);
        rightDists = grid->globalGridDims - grid->numInternalCells;
        n = 4;
        const FP r0 = 1e-8;
        for (int d = 0; d < 3; d++)
            if (globalLeftDims[d])
                maxSigma[d] = -log(r0) * (n + 1) * constants::c / // constants::c is necessary because different cgs are used
                (2 * globalLeftDims[d] * grid->steps[d]);
            else
                maxSigma[d] = 0;
    }

    template<GridTypes gridTypes>
    FP3 Pml<gridTypes>::computeSigma(const FP3& coords)
    {
        Grid<FP, gridTypes>* grid = fieldSolver->grid;
        FP3 globalLeftBorder = grid->origin - grid->steps * (leftDists);
        FP3 pmlLeftBorder = globalLeftBorder + grid->steps * (globalLeftDims + grid->getNumExternalLeftCells());
        FP3 globalRightBorder = globalLeftBorder + grid->steps *
            (leftDists + rightDists + grid->numCells);
        FP3 pmlRightBorder = globalRightBorder - grid->steps * (globalRightDims + grid->getNumExternalRightCells());

        // for each dimension coeff is coefficient of linear interpolation between 0
        // and 1 with 0 on PML-internal area boundary and 1 on outer PML boundary,
        // sigma = maxSigma * coeff^n
        FP3 sigma;
        for (int d = 0; d < 3; ++d)
            if (coords[d] < pmlLeftBorder[d])
            {
                FP coeff = (pmlLeftBorder[d] - coords[d]) /
                    (pmlLeftBorder[d] - globalLeftBorder[d]);
                sigma[d] = maxSigma[d] * pow(coeff, n);
            }
            else
                if (coords[d] > pmlRightBorder[d])
                {
                    FP coeff = (coords[d] - pmlRightBorder[d]) /
                        (globalRightBorder[d] - pmlRightBorder[d]);
                    sigma[d] = maxSigma[d] * pow(coeff, n);
                }
        // else sigma[d] = 0 as it is already
        return sigma;
    }

    template<GridTypes gridTypes>
    void Pml<gridTypes>::initializeSplitFieldsE(int sizeE) {
        exy.resize(sizeE, 0);
        exz.resize(sizeE, 0);
        eyx.resize(sizeE, 0);
        eyz.resize(sizeE, 0);
        ezx.resize(sizeE, 0);
        ezy.resize(sizeE, 0);
    }

    template<GridTypes gridTypes>
    void Pml<gridTypes>::initializeSplitFieldsB(int sizeB) {
        bxy.resize(sizeB, 0);
        bxz.resize(sizeB, 0);
        byx.resize(sizeB, 0);
        byz.resize(sizeB, 0);
        bzx.resize(sizeB, 0);
        bzy.resize(sizeB, 0);
    }


    template<GridTypes gridTypes>
    class PmlReal : public Pml<gridTypes>
    {
    public:

        PmlReal(RealFieldSolver<gridTypes>* fieldSolver, Int3 sizePML);

        void computeCoeffs();
        void updateDims();

        virtual void updateB() {};
        virtual void updateE() {};

        int numNodes, numCells; // total number of PML nodes / cells
        std::vector<Int3> nodeIndex, cellIndex; // natural 3d indexes of nodes in PML
        std::vector<int> spaceIndexToPmlCell, spaceIndexToPmlNode;
        std::vector<FP3> coeffEa, coeffEb, coeffBa, coeffBb; // coeffs for FDTD in PML
        std::vector<FP> coeffJ;
    };

    template<GridTypes gridTypes>
    PmlReal<gridTypes>::PmlReal(RealFieldSolver<gridTypes>* _fieldSolver, Int3 _sizePML) :
        Pml<gridTypes>((FieldSolver<gridTypes>*) _fieldSolver, _sizePML)
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;

        int totalStored = grid->numCells.x * grid->numCells.y * grid->numCells.z;
        spaceIndexToPmlNode = std::vector<int>(totalStored, -1);
        spaceIndexToPmlCell = std::vector<int>(totalStored, -1);

        // Fill data of PML values of electric field
        const Int3 leftPmlEnd = this->leftDims;
        const Int3 rightPmlBegin = grid->numCells - this->rightDims;
        const int maxSize = grid->numCells.volume() -
            (grid->numCells - this->leftDims - this->rightDims).volume();
        Int3 begin = this->fieldSolver->updateEAreaBegin;
        Int3 end = this->fieldSolver->updateEAreaEnd + Int3(1, 1, 1); // + 1 for edge values of E
        for (int i = grid->dimensionality; i < 3; i++)
            end[i]--;
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
                    bool xBoundaryPml = (i < leftPmlEnd.x) || (i >= rightPmlBegin.x);
                    bool yBoundaryPml = (j < leftPmlEnd.y) || (j >= rightPmlBegin.y);
                    bool zBoundaryPml = (k < leftPmlEnd.z) || (k >= rightPmlBegin.z);
                    if (xBoundaryPml || yBoundaryPml || zBoundaryPml)
                    {
                        nodeIndex.push_back(Int3(i, j, k));
                        int index = k + (j + i * grid->numCells.y) * grid->numCells.z;
                        spaceIndexToPmlNode[index] = (int)nodeIndex.size() - 1;
                    }
                }

        numNodes = (int)nodeIndex.size();
        this->initializeSplitFieldsE(numNodes);

        // Fill data of PML values of magnetic field
        begin = this->fieldSolver->updateBAreaBegin;
        end = this->fieldSolver->updateBAreaEnd;
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
                        spaceIndexToPmlCell[index] = (int)cellIndex.size() - 1;
                    }
                }
        numCells = (int)cellIndex.size();
        this->initializeSplitFieldsB(numCells);

        coeffEa.resize(numNodes);
        coeffEb.resize(numNodes);
        coeffBa.resize(numCells);
        coeffBb.resize(numCells);
        coeffJ.resize(numNodes);
        computeCoeffs();
    }

    template<GridTypes gridTypes>
    void PmlReal<gridTypes>::computeCoeffs()
    {
        Grid<FP, gridTypes> * grid = this->fieldSolver->grid;
        FP cdt = constants::c * this->fieldSolver->dt;
        const FP threshold = (FP)1e-8;
        coeffEa.resize(numNodes);
        coeffEb.resize(numNodes);
        coeffBa.resize(numCells);
        coeffBb.resize(numCells);
        coeffJ.resize(numNodes);

        for (int idx = 0; idx < numNodes; ++idx)
        {
            int i = nodeIndex[idx].x;
            int j = nodeIndex[idx].y;
            int k = nodeIndex[idx].z;
            FP3 eCoords[] = { grid->EyPosition(i, j, k), grid->ExPosition(i, j, k),
                grid->ExPosition(i, j, k) };
            for (int d = 0; d < 3; d++)
            {
                FP3 sigma = this->computeSigma(eCoords[d]);
                int axis1 = (d + 1) % 3, axis2 = (d + 2) % 3;
                if (d == 1)
                {
                    axis1 = 0;
                    axis2 = 2;
                }
                if (sigma[d] > threshold)
                {
                    FP expCoeff = exp(-sigma[d] * cdt);
                    coeffEa[idx][d] = expCoeff;
                    coeffEb[idx][d] = (expCoeff - 1) / (sigma[d] * grid->steps[d]);
                }
                else
                {
                    coeffEa[idx][d] = 1;
                    coeffEb[idx][d] = -cdt / (grid->steps[d]);
                }
            }
            coeffJ[idx] = (FP)1;
            for (int d = 0; d < 3; d++)
                if (coeffEa[idx][d] < coeffJ[idx])
                    coeffJ[idx] = coeffEa[idx][d];
            coeffJ[idx] *= -(FP)4 * constants::pi * this->fieldSolver->dt / (FP)2;
            // divide over 2 because of half fields
        }

        for (int idx = 0; idx < numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;
            FP3 bCoords[] = { grid->ByPosition(i, j, k), grid->BxPosition(i, j, k),
                grid->BxPosition(i, j, k) };
            for (int d = 0; d < 3; d++)
            {
                FP3 sigma = this->computeSigma(bCoords[d]);
                int axis1 = (d + 1) % 3, axis2 = (d + 2) % 3;
                if (d == 1)
                {
                    axis1 = 0;
                    axis2 = 2;
                }
                if (sigma[d] > threshold)
                {
                    FP expCoeff = exp(-sigma[d] * cdt);
                    coeffBa[idx][d] = expCoeff;
                    coeffBb[idx][d] = -(expCoeff - 1) / (sigma[d] * grid->steps[d]);
                }
                else
                {
                    coeffBa[idx][d] = 1;
                    coeffBb[idx][d] = cdt / (grid->steps[d]);
                }
            }
        }
    }

    template<GridTypes gridTypes>
    class PmlSpectral : public Pml<gridTypes>
    {
    public:

        PmlSpectral(SpectralFieldSolver<gridTypes>* fieldSolver, Int3 sizePML);

        void computeCoeffs();

        // first step of spectral pml: updating of split components
        virtual void updateExSplit() {};
        virtual void updateEySplit() {};
        virtual void updateEzSplit() {};
        virtual void updateBxSplit() {};
        virtual void updateBySplit() {};
        virtual void updateBzSplit() {};

        virtual void updateBSplit();
        virtual void updateESplit();

        // first step of spectral pml: updating of split components
        virtual void doFirstStep();

        // multiplication by exp(-sigma*dt)
        void multBySigmaE();
        void multBySigmaB();
        void multBySigma();

        // summation of split components
        void sumField(ScalarField<FP>& field, std::vector<FP>& splitField1,
            std::vector<FP>& splitField2);
        void sumFieldE();
        void sumFieldB();
        void sumFields();

        // second step of spectral pml: multiplication by exp(-sigma*dt) and summation of split components
        void doSecondStep();
        
        int numCells; // total number of PML cells
        std::vector<Int3> cellIndex; // natural 3d indexes of cells in PML
        std::vector<int> spaceIndexToPmlCell;
        std::vector<FP3> coeff;

    };

    template<GridTypes gridTypes>
    inline PmlSpectral<gridTypes>::PmlSpectral(SpectralFieldSolver<gridTypes>* _fieldSolver, Int3 _sizePML) :
        Pml<gridTypes>((FieldSolver<gridTypes>*) _fieldSolver, _sizePML)
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;

        spaceIndexToPmlCell = std::vector<int>(grid->numCells.volume(), -1);

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
                        spaceIndexToPmlCell[index] = (int)cellIndex.size() - 1;
                    }
                }

        numCells = (int)cellIndex.size();
        this->initializeSplitFieldsE(numCells);
        this->initializeSplitFieldsB(numCells);

        computeCoeffs();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::computeCoeffs()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        coeff.resize(numCells);
        for (int idx = 0; idx < numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            FP3 coords = grid->ExPosition(i, j, k);
            FP3 sigma = this->computeSigma(coords);

            coeff[idx].x = exp(-sigma.x * this->fieldSolver->dt);
            coeff[idx].y = exp(-sigma.y * this->fieldSolver->dt);
            coeff[idx].z = exp(-sigma.z * this->fieldSolver->dt);
        }
    }
    
    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::updateBSplit()
    {
        updateBxSplit();
        updateBySplit();
        updateBzSplit();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::updateESplit()
    {
        updateExSplit();
        updateEySplit();
        updateEzSplit();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::doFirstStep()
    {
        updateBSplit();
        updateESplit();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::multBySigmaE()
    {
        for (int idx = 0; idx < numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            this->eyx[idx] *= coeff[idx].y;
            this->ezx[idx] *= coeff[idx].z;
            this->ezy[idx] *= coeff[idx].z;
            this->exy[idx] *= coeff[idx].x;
            this->exz[idx] *= coeff[idx].x;
            this->eyz[idx] *= coeff[idx].y;
        }
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::multBySigmaB()
    {
        for (int idx = 0; idx < numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            this->byx[idx] *= coeff[idx].y;
            this->bzx[idx] *= coeff[idx].z;
            this->bzy[idx] *= coeff[idx].z;
            this->bxy[idx] *= coeff[idx].x;
            this->bxz[idx] *= coeff[idx].x;
            this->byz[idx] *= coeff[idx].y;
        }
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::multBySigma()
    {
        multBySigmaE();
        multBySigmaB();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumField(ScalarField<FP>& field, std::vector<FP>& splitField1,
        std::vector<FP>& splitField2)
    {
        for (int idx = 0; idx < numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

           field(i, j, k) = splitField1[idx] + splitField2[idx];
        }
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldE()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        sumField(grid->Ex, this->eyx, this->ezx);
        sumField(grid->Ey, this->ezy, this->exy);
        sumField(grid->Ez, this->exz, this->eyz);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldB()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        sumField(grid->Bx, this->byx, this->bzx);
        sumField(grid->By, this->bzy, this->bxy);
        sumField(grid->Bz, this->bxz, this->byz);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFields()
    {
        sumFieldE();
        sumFieldB();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::doSecondStep()
    {
        multBySigma();
        sumFields();
    }
}