#pragma once
#include <vector>

#include "Grid.h"
#include "PmlSplitGrid.h"

namespace pfc {

    template<class TGrid>
    class Pml
    {
    public:

        Pml(TGrid* grid, FP dt, Int3 domainIndexBegin, Int3 domainIndexEnd,
            Int3 sizePML, FP nPmlParam, FP r0PmlParam);

        // constructor for loading
        explicit Pml(TGrid* grid, FP dt, Int3 domainIndexBegin, Int3 domainIndexEnd);

        /* implement the next methods in derived classes
        void updateB();
        void updateE();
        */

        forceinline FP computeSigma(FP coord, CoordinateEnum axis) const;

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        TGrid* grid = nullptr;
        FP dt = 0.0;
        Int3 domainIndexBegin, domainIndexEnd;  // internal domain area including pml

        std::unique_ptr<PmlSplitGrid> splitGrid;

        Int3 sizePML;
        Int3 leftPmlBorder, rightPmlBorder;
        FP3 leftPmlBorderCoord, rightPmlBorderCoord;
        FP3 leftGlobalBorderCoord, rightGlobalBorderCoord;

    protected:

        // absorbing function parameters
        FP3 maxSigma;
        FP n = 0.0;

        void initializePmlParams(FP n, FP r0, const FP3& gridStep);

        void checkGridAndPmlSize(const Int3& globalGridDims) {
            if (2 * sizePML > globalGridDims) {
                std::string exc = "ERROR: this->grid size should be larger than 2x pml size";
                std::cout << exc << std::endl;
                throw std::logic_error(exc);
            }
        }

    };

    template<class TGrid>
    inline void Pml<TGrid>::initializePmlParams(FP n, FP r0, const FP3& gridStep) {
        this->n = n;

        for (int d = 0; d < 3; d++) {
            if (this->sizePML[d])
                this->maxSigma[d] = -log(r0) * (n + 1) / (2 * this->sizePML[d] * gridStep[d]);
            else this->maxSigma[d] = 0;
        }
    }

    template<class TGrid>
    inline Pml<TGrid>::Pml(TGrid* grid, FP dt,
        Int3 domainIndexBegin, Int3 domainIndexEnd) :
        grid(grid), dt(dt),
        domainIndexBegin(domainIndexBegin), domainIndexEnd(domainIndexEnd)
    {}

    template<class TGrid>
    inline Pml<TGrid>::Pml(TGrid* grid, FP dt,
        Int3 domainIndexBegin, Int3 domainIndexEnd,
        Int3 sizePML, FP nPmlParam, FP r0PmlParam) :
        grid(grid), dt(dt), sizePML(sizePML),
        domainIndexBegin(domainIndexBegin), domainIndexEnd(domainIndexEnd)
    {
        checkGridAndPmlSize(this->grid->globalGridDims);

        // TODO: consider global and local domain borders
        this->leftPmlBorder = this->sizePML + this->domainIndexBegin;
        this->rightPmlBorder = this->domainIndexEnd - this->sizePML;
        this->leftGlobalBorderCoord = this->grid->origin + (FP3)this->domainIndexBegin * this->grid->steps;
        this->rightGlobalBorderCoord = this->grid->origin + (FP3)this->domainIndexEnd * this->grid->steps;
        this->leftPmlBorderCoord = this->leftGlobalBorderCoord + this->grid->steps * (FP3)this->sizePML;
        this->rightPmlBorderCoord = this->rightGlobalBorderCoord - this->grid->steps * (FP3)this->sizePML;

        this->splitGrid.reset(new PmlSplitGrid(this->leftPmlBorder, this->rightPmlBorder,
            this->domainIndexBegin, this->domainIndexEnd));
        initializePmlParams(nPmlParam, r0PmlParam, this->grid->steps);
    }

    template<class TGrid>
    inline FP Pml<TGrid>::computeSigma(FP coord, CoordinateEnum axis) const
    {
        int d = (int)axis;

        FP distance;
        if (coord < leftPmlBorderCoord[d])
            distance = (leftPmlBorderCoord[d] - coord) / (leftPmlBorderCoord[d] - leftGlobalBorderCoord[d]);
        else if (coord >= rightPmlBorderCoord[d])
            distance = (coord - rightPmlBorderCoord[d]) / (rightGlobalBorderCoord[d] - rightPmlBorderCoord[d]);
        else distance = 0;

        return this->maxSigma[d] * pow(distance, (FP)n);
    }

    template<class TGrid>
    inline void Pml<TGrid>::save(std::ostream& ostr)
    {
        ostr.write((char*)&sizePML, sizeof(sizePML));
        ostr.write((char*)&n, sizeof(n));
        ostr.write((char*)&maxSigma, sizeof(maxSigma));
        ostr.write((char*)&leftPmlBorder, sizeof(leftPmlBorder));
        ostr.write((char*)&rightPmlBorder, sizeof(rightPmlBorder));
        ostr.write((char*)&leftPmlBorderCoord, sizeof(leftPmlBorderCoord));
        ostr.write((char*)&rightPmlBorderCoord, sizeof(rightPmlBorderCoord));
        ostr.write((char*)&leftGlobalBorderCoord, sizeof(leftGlobalBorderCoord));
        ostr.write((char*)&rightGlobalBorderCoord, sizeof(rightGlobalBorderCoord));

        splitGrid->save(ostr);
    }

    template<class TGrid>
    inline void Pml<TGrid>::load(std::istream& istr)
    {
        istr.read((char*)&sizePML, sizeof(sizePML));
        istr.read((char*)&n, sizeof(n));
        istr.read((char*)&maxSigma, sizeof(maxSigma));
        istr.read((char*)&leftPmlBorder, sizeof(leftPmlBorder));
        istr.read((char*)&rightPmlBorder, sizeof(rightPmlBorder));
        istr.read((char*)&leftPmlBorderCoord, sizeof(leftPmlBorderCoord));
        istr.read((char*)&rightPmlBorderCoord, sizeof(rightPmlBorderCoord));
        istr.read((char*)&leftGlobalBorderCoord, sizeof(leftGlobalBorderCoord));
        istr.read((char*)&rightGlobalBorderCoord, sizeof(rightGlobalBorderCoord));

        this->splitGrid.reset(new PmlSplitGrid(this->leftPmlBorder, this->rightPmlBorder,
            this->domainIndexBegin, this->domainIndexEnd));
        splitGrid->load(istr);
    }


    template<class TGrid>
    class PmlReal : public Pml<TGrid>
    {
    public:

        PmlReal(TGrid* grid, FP dt, Int3 domainIndexBegin, Int3 domainIndexEnd,
            Int3 sizePML, FP nPmlParam, FP r0PmlParam);

        // constructor for loading
        explicit PmlReal(TGrid* grid, FP dt, Int3 domainIndexBegin, Int3 domainIndexEnd);

    };

    template<class TGrid>
    inline PmlReal<TGrid>::PmlReal(TGrid* grid, FP dt, Int3 domainIndexBegin, Int3 domainIndexEnd,
        Int3 sizePML, FP nPmlParam, FP r0PmlParam) :
        Pml<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd,
            sizePML, nPmlParam, r0PmlParam)
    {}

    template<class TGrid>
    inline PmlReal<TGrid>::PmlReal(TGrid* grid, FP dt, Int3 domainIndexBegin, Int3 domainIndexEnd) :
        Pml<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd)
    {}


    template<class TGrid>
    class PmlSpectral : public Pml<TGrid>
    {
    public:

        PmlSpectral(TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt,
            Int3 domainIndexBegin, Int3 domainIndexEnd,
            Int3 complexDomainIndexBegin, Int3 complexDomainIndexEnd,
            Int3 sizePML, FP nPmlParam, FP r0PmlParam);

        // constructor for loading
        explicit PmlSpectral(TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt,
            Int3 domainIndexBegin, Int3 domainIndexEnd,
            Int3 complexDomainIndexBegin, Int3 complexDomainIndexEnd);

        /* implement the next methods in derived classes
        // first step of spectral pml: updating of split components
        void updateBSplit();
        void updateESplit();
        */

        // second step of spectral pml: multiplication by exp(-sigma*dt) and summation of split components
        void updateB();
        void updateE();

        SpectralGrid<FP, complexFP>* complexGrid = nullptr;
        Int3 complexDomainIndexBegin, complexDomainIndexEnd;
    };

    template<class TGrid>
    inline PmlSpectral<TGrid>::PmlSpectral(
        TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt,
        Int3 domainIndexBegin, Int3 domainIndexEnd,
        Int3 complexDomainIndexBegin, Int3 complexDomainIndexEnd,
        Int3 sizePML, FP nPmlParam, FP r0PmlParam) :
        Pml<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd,
            sizePML, nPmlParam, r0PmlParam),
        complexGrid(complexGrid),
        complexDomainIndexBegin(complexDomainIndexBegin),
        complexDomainIndexEnd(complexDomainIndexEnd)
    {}

    template<class TGrid>
    inline PmlSpectral<TGrid>::PmlSpectral(
        TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt,
        Int3 domainIndexBegin, Int3 domainIndexEnd,
        Int3 complexDomainIndexBegin, Int3 complexDomainIndexEnd) :
        Pml<TGrid>(grid, dt, domainIndexBegin, domainIndexEnd),
        complexGrid(complexGrid),
        complexDomainIndexBegin(complexDomainIndexBegin),
        complexDomainIndexEnd(complexDomainIndexEnd)
    {}

    template<class TGrid>
    inline void PmlSpectral<TGrid>::updateB()
    {
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index = this->splitGrid->getIndex3d(idx);

            FP sigmaX = this->computeSigma(this->grid->ByPosition(index.x, index.y, index.z).x, CoordinateEnum::x);
            FP sigmaY = this->computeSigma(this->grid->BzPosition(index.x, index.y, index.z).y, CoordinateEnum::y);
            FP sigmaZ = this->computeSigma(this->grid->BxPosition(index.x, index.y, index.z).z, CoordinateEnum::z);
            FP coeffX = exp(-sigmaX * cdt), coeffY = exp(-sigmaY * cdt), coeffZ = exp(-sigmaZ * cdt);

            this->splitGrid->bxy[idx] *= coeffY;
            this->splitGrid->bxz[idx] *= coeffZ;
            this->splitGrid->byz[idx] *= coeffZ;
            this->splitGrid->byx[idx] *= coeffX;
            this->splitGrid->bzx[idx] *= coeffX;
            this->splitGrid->bzy[idx] *= coeffY;

            this->grid->Bx(index.x, index.y, index.z) = this->splitGrid->bxy[idx] + this->splitGrid->bxz[idx];
            this->grid->By(index.x, index.y, index.z) = this->splitGrid->byz[idx] + this->splitGrid->byx[idx];
            this->grid->Bz(index.x, index.y, index.z) = this->splitGrid->bzx[idx] + this->splitGrid->bzy[idx];
        }
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::updateE()
    {
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index = this->splitGrid->getIndex3d(idx);

            FP sigmaX = this->computeSigma(this->grid->EyPosition(index.x, index.y, index.z).x, CoordinateEnum::x);
            FP sigmaY = this->computeSigma(this->grid->EzPosition(index.x, index.y, index.z).y, CoordinateEnum::y);
            FP sigmaZ = this->computeSigma(this->grid->ExPosition(index.x, index.y, index.z).z, CoordinateEnum::z);
            FP coeffX = exp(-sigmaX * cdt), coeffY = exp(-sigmaY * cdt), coeffZ = exp(-sigmaZ * cdt);

            this->splitGrid->exy[idx] *= coeffY;
            this->splitGrid->exz[idx] *= coeffZ;
            this->splitGrid->eyz[idx] *= coeffZ;
            this->splitGrid->eyx[idx] *= coeffX;
            this->splitGrid->ezx[idx] *= coeffX;
            this->splitGrid->ezy[idx] *= coeffY;

            this->grid->Ex(index.x, index.y, index.z) = this->splitGrid->exy[idx] + this->splitGrid->exz[idx];
            this->grid->Ey(index.x, index.y, index.z) = this->splitGrid->eyz[idx] + this->splitGrid->eyx[idx];
            this->grid->Ez(index.x, index.y, index.z) = this->splitGrid->ezx[idx] + this->splitGrid->ezy[idx];
        }
    }
}
