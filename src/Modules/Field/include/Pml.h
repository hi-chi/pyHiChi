#pragma once
#include <vector>

#include "Grid.h"

namespace pfc {

    template<class TGrid>
    class Pml
    {
    public:
        Pml(TGrid* grid, FP dt, Int3 sizePML,
            FP nPmlParam = (FP)4.0, FP r0PmlParam = (FP)1e-8);

        /* implement the next methods in derived classes
        void updateB();
        void updateE();
        */

        forceinline FP computeSigma(FP coord, CoordinateEnum axis) const;

        TGrid* grid = nullptr;
        FP dt = 0.0;

        Int3 sizePML;
        Int3 leftPmlBorder, rightPmlBorder;
        std::vector<FP> bxy, bxz, byx, byz, bzx, bzy;  // split magnetic field
        std::vector<FP> exy, exz, eyx, eyz, ezx, ezy;  // split electric field

        std::vector<Int3> bIndex, eIndex;  // natural 3d indexes of nodes in PML

    protected:

        // absorbing function parameters
        FP3 maxSigma;
        FP n;

        void initializePmlParams(FP _n, FP _r0, const FP3& gridStep);

        void initializeSplitFieldsE(int sizeE);
        void initializeSplitFieldsB(int sizeB);

        void fillPmlIndices();

        void checkGridAndPmlSize(const Int3& globalGridDims) {
            if (2 * sizePML > globalGridDims) {
                std::string exc = "ERROR: this->grid size should be larger than 2x pml size";
                std::cout << exc << std::endl;
                throw std::logic_error(exc);
            }
        }

    };

    template<class TGrid>
    void Pml<TGrid>::initializePmlParams(FP n, FP r0, const FP3& gridStep) {
        this->n = n;

        for (int d = 0; d < 3; d++) {
            if (this->sizePML[d])
                this->maxSigma[d] = -log(r0) * (n + 1) / (2 * this->sizePML[d] * gridStep[d]);
            else this->maxSigma[d] = 0;
        }
    }

    template<class TGrid>
    inline Pml<TGrid>::Pml(TGrid* grid, FP dt, Int3 sizePML,
        FP nPmlParam, FP r0PmlParam) :
        grid(grid), dt(dt), sizePML(sizePML)
    {
        checkGridAndPmlSize(this->grid->globalGridDims);

        // TODO: consider that pml is only in edge domains
        leftPmlBorder = this->sizePML + this->grid->getNumExternalLeftCells();
        rightPmlBorder = leftPmlBorder + (this->grid->numInternalCells - 2 * this->sizePML);

        initializePmlParams(nPmlParam, r0PmlParam, this->grid->steps);
        fillPmlIndices();
    }

    template<class TGrid>
    void Pml<TGrid>::fillPmlIndices()
    {
        // TODO: check border indices
        Int3 begin = this->grid->correctNumCellsAccordingToDim(Int3(1, 1, 1));
        Int3 end = this->grid->numCells;

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
                    bool xBoundaryPml = (i < this->leftPmlBorder.x) || (i >= this->rightPmlBorder.x);
                    bool yBoundaryPml = (j < this->leftPmlBorder.y) || (j >= this->rightPmlBorder.y);
                    bool zBoundaryPml = (k < this->leftPmlBorder.z) || (k >= this->rightPmlBorder.z);
                    if (xBoundaryPml || yBoundaryPml || zBoundaryPml)
                    {
                        bIndex.push_back(Int3(i, j, k));
                    }
                }

        begin = Int3(0, 0, 0);
        end = this->grid->numCells -
            this->grid->correctNumCellsAccordingToDim(Int3(1, 1, 1));;

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
                    bool xBoundaryPml = (i < this->leftPmlBorder.x) || (i >= this->rightPmlBorder.x);
                    bool yBoundaryPml = (j < this->leftPmlBorder.y) || (j >= this->rightPmlBorder.y);
                    bool zBoundaryPml = (k < this->leftPmlBorder.z) || (k >= this->rightPmlBorder.z);
                    if (xBoundaryPml || yBoundaryPml || zBoundaryPml)
                    {
                        eIndex.push_back(Int3(i, j, k));
                    }
                }

        this->initializeSplitFieldsB(bIndex.size());
        this->initializeSplitFieldsE(eIndex.size());
    }

    template<class TGrid>
    void Pml<TGrid>::initializeSplitFieldsB(int sizeB) {
        bxy.resize(sizeB, 0);
        bxz.resize(sizeB, 0);
        byx.resize(sizeB, 0);
        byz.resize(sizeB, 0);
        bzx.resize(sizeB, 0);
        bzy.resize(sizeB, 0);
    }

    template<class TGrid>
    void Pml<TGrid>::initializeSplitFieldsE(int sizeE) {
        exy.resize(sizeE, 0);
        exz.resize(sizeE, 0);
        eyx.resize(sizeE, 0);
        eyz.resize(sizeE, 0);
        ezx.resize(sizeE, 0);
        ezy.resize(sizeE, 0);
    }


    template<class TGrid>
    FP Pml<TGrid>::computeSigma(FP coord, CoordinateEnum axis) const
    {
        int d = (int)axis;

        FP globalLeftBorderCoord = this->grid->origin[d];
        FP globalRightBorderCoord = this->grid->origin[d] + this->grid->steps[d] * this->grid->numCells[d];
        FP leftPmlBorderCoord = this->grid->origin[d] + this->grid->steps[d] * leftPmlBorder[d];
        FP rightPmlBorderCoord = this->grid->origin[d] + this->grid->steps[d] * rightPmlBorder[d];

        FP distance;
        if (coord < leftPmlBorderCoord)
            distance = (leftPmlBorderCoord - coord) / (leftPmlBorderCoord - globalLeftBorderCoord);
        else if (coord >= rightPmlBorderCoord)
            distance = (coord - rightPmlBorderCoord) / (globalRightBorderCoord - rightPmlBorderCoord);
        else distance = 0;

        return this->maxSigma[d] * pow(distance, (FP)n);
    }


    template<class TGrid>
    class PmlReal : public Pml<TGrid>
    {
    public:

        PmlReal(TGrid* grid, FP dt, Int3 sizePML,
            FP nPmlParam = (FP)4.0, FP r0PmlParam = (FP)1e-8) :
            Pml<TGrid>(grid, dt, sizePML, nPmlParam, r0PmlParam)
        {}

    };


    template<class TGrid>
    class PmlSpectral : public Pml<TGrid>
    {
    public:

        PmlSpectral(TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt,
            Int3 sizePML, FP nPmlParam = (FP)4.0, FP r0PmlParam = (FP)1e-8) :
            Pml<TGrid>(grid, dt, sizePML, nPmlParam, r0PmlParam),
            complexGrid(complexGrid)
        {}

        /* implement the next methods in derived classes
        // first step of spectral pml: updating of split components
        void updateExSplit();
        void updateEySplit();
        void updateEzSplit();
        void updateBxSplit();
        void updateBySplit();
        void updateBzSplit();

        void updateBSplit();
        void updateESplit();
        */

        // multiplication by exp(-sigma*dt)
        void multBySigmaE();
        void multBySigmaB();

        // summation of split components
        void sumFieldEx();
        void sumFieldEy();
        void sumFieldEz();
        void sumFieldBx();
        void sumFieldBy();
        void sumFieldBz();

        void sumFieldE();
        void sumFieldB();

        // second step of spectral pml: multiplication by exp(-sigma*dt) and summation of split components
        void updateB();
        void updateE();

        SpectralGrid<FP, complexFP>* complexGrid = nullptr;

    protected:

        void sumFieldComponent(ScalarField<FP>& field, const std::vector<Int3>& index,
            std::vector<FP>& splitField1, std::vector<FP>& splitField2);
    };

    template<class TGrid>
    inline void PmlSpectral<TGrid>::multBySigmaE()
    {
        const int size = this->eIndex.size();
        const FP cdt = constants::c * this->dt;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            int i = this->eIndex[idx].x, j = this->eIndex[idx].y, k = this->eIndex[idx].z;

            FP sigmaX = this->computeSigma(this->grid->EyPosition(i, j, k).x, CoordinateEnum::x);
            FP sigmaY = this->computeSigma(this->grid->EzPosition(i, j, k).y, CoordinateEnum::y);
            FP sigmaZ = this->computeSigma(this->grid->ExPosition(i, j, k).z, CoordinateEnum::z);
            FP coeffX = exp(-sigmaX * cdt), coeffY = exp(-sigmaY * cdt), coeffZ = exp(-sigmaZ * cdt);

            this->eyx[idx] *= coeffY;
            this->ezx[idx] *= coeffZ;
            this->ezy[idx] *= coeffZ;
            this->exy[idx] *= coeffX;
            this->exz[idx] *= coeffX;
            this->eyz[idx] *= coeffY;
        }
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::multBySigmaB()
    {
        const int size = this->bIndex.size();
        const FP cdt = constants::c * this->dt;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            int i = this->bIndex[idx].x, j = this->bIndex[idx].y, k = this->bIndex[idx].z;

            FP sigmaX = this->computeSigma(this->grid->ByPosition(i, j, k).x, CoordinateEnum::x);
            FP sigmaY = this->computeSigma(this->grid->BzPosition(i, j, k).y, CoordinateEnum::y);
            FP sigmaZ = this->computeSigma(this->grid->BxPosition(i, j, k).z, CoordinateEnum::z);
            FP coeffX = exp(-sigmaX * cdt), coeffY = exp(-sigmaY * cdt), coeffZ = exp(-sigmaZ * cdt);

            this->byx[idx] *= coeffY;
            this->bzx[idx] *= coeffZ;
            this->bzy[idx] *= coeffZ;
            this->bxy[idx] *= coeffX;
            this->bxz[idx] *= coeffX;
            this->byz[idx] *= coeffY;
        }
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::sumFieldEx()
    {
        sumFieldComponent(this->grid->Ex, this->eIndex, this->eyx, this->ezx);
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::sumFieldEy()
    {
        sumFieldComponent(this->grid->Ey, this->eIndex, this->ezy, this->exy);
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::sumFieldEz()
    {
        sumFieldComponent(this->grid->Ez, this->eIndex, this->exz, this->eyz);
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::sumFieldE()
    {
        this->sumFieldEx();
        this->sumFieldEy();
        this->sumFieldEz();
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::sumFieldBx()
    {
        sumFieldComponent(this->grid->Bx, this->bIndex, this->byx, this->bzx);
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::sumFieldBy()
    {
        sumFieldComponent(this->grid->By, this->bIndex, this->bzy, this->bxy);
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::sumFieldBz()
    {
        sumFieldComponent(this->grid->Bz, this->bIndex, this->bxz, this->byz);
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::sumFieldB()
    {
        this->sumFieldBx();
        this->sumFieldBy();
        this->sumFieldBz();
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::updateB()
    {
        multBySigmaB();
        sumFieldB();
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::updateE()
    {
        multBySigmaE();
        sumFieldE();
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::sumFieldComponent(ScalarField<FP>& field,
        const std::vector<Int3>& index,
        std::vector<FP>& splitField1, std::vector<FP>& splitField2)
    {
        const int size = index.size();

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            int i = index[idx].x, j = index[idx].y, k = index[idx].z;

            field(i, j, k) = splitField1[idx] + splitField2[idx];
        }
    }
}
