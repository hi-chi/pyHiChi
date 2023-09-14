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

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        // coefficient pre-computing
        void computeCoeffs();

        SpectralGrid<FP, complexFP>* complexGrid = nullptr;
        Int3 complexDomainIndexBegin, complexDomainIndexEnd;

        std::vector<FP> bCoeffX, bCoeffY, bCoeffZ, eCoeffX, eCoeffY, eCoeffZ;

    private:

        void computeCoeffs(
            std::vector<FP>& coeffX, std::vector<FP>& coeffY, std::vector<FP>& coeffZ,
            const FP3(TGrid::* positionFX)(int, int, int) const,
            const FP3(TGrid::* positionFY)(int, int, int) const,
            const FP3(TGrid::* positionFZ)(int, int, int) const);

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
    {
        this->computeCoeffs();
    }

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
    inline void PmlSpectral<TGrid>::computeCoeffs()
    {
        this->computeCoeffs(bCoeffX, bCoeffY, bCoeffZ,
            &TGrid::BxPosition, &TGrid::ByPosition, &TGrid::BzPosition);
        this->computeCoeffs(eCoeffX, eCoeffY, eCoeffZ,
            &TGrid::ExPosition, &TGrid::EyPosition, &TGrid::EzPosition);
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::updateB()
    {
        const int size = this->splitGrid->getNumPmlNodes();

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index = this->splitGrid->getIndex3d(idx);

            this->splitGrid->bxy[idx] *= bCoeffY[idx];
            this->splitGrid->bxz[idx] *= bCoeffZ[idx];
            this->splitGrid->byz[idx] *= bCoeffZ[idx];
            this->splitGrid->byx[idx] *= bCoeffX[idx];
            this->splitGrid->bzx[idx] *= bCoeffX[idx];
            this->splitGrid->bzy[idx] *= bCoeffY[idx];

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

            this->splitGrid->exy[idx] *= eCoeffY[idx];
            this->splitGrid->exz[idx] *= eCoeffZ[idx];
            this->splitGrid->eyz[idx] *= eCoeffZ[idx];
            this->splitGrid->eyx[idx] *= eCoeffX[idx];
            this->splitGrid->ezx[idx] *= eCoeffX[idx];
            this->splitGrid->ezy[idx] *= eCoeffY[idx];

            this->grid->Ex(index.x, index.y, index.z) = this->splitGrid->exy[idx] + this->splitGrid->exz[idx];
            this->grid->Ey(index.x, index.y, index.z) = this->splitGrid->eyz[idx] + this->splitGrid->eyx[idx];
            this->grid->Ez(index.x, index.y, index.z) = this->splitGrid->ezx[idx] + this->splitGrid->ezy[idx];
        }
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::computeCoeffs(
        std::vector<FP>& coeffX, std::vector<FP>& coeffY, std::vector<FP>& coeffZ,
        const FP3(TGrid::* positionFX)(int, int, int) const,
        const FP3(TGrid::* positionFY)(int, int, int) const,
        const FP3(TGrid::* positionFZ)(int, int, int) const)
    {
        const int size = this->splitGrid->getNumPmlNodes();

        coeffX.resize(size);
        coeffY.resize(size);
        coeffZ.resize(size);

        const FP cdt = constants::c * this->dt;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index = this->splitGrid->getIndex3d(idx);

            // coordinates according to Yee grid or collocated grid
            FP sigmaX = this->computeSigma((this->grid->*positionFY)(index.x, index.y, index.z).x, CoordinateEnum::x);
            FP sigmaY = this->computeSigma((this->grid->*positionFZ)(index.x, index.y, index.z).y, CoordinateEnum::y);
            FP sigmaZ = this->computeSigma((this->grid->*positionFX)(index.x, index.y, index.z).z, CoordinateEnum::z);

            coeffX[idx] = exp(-sigmaX * cdt);
            coeffY[idx] = exp(-sigmaY * cdt);
            coeffZ[idx] = exp(-sigmaZ * cdt);
        }
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::save(std::ostream& ostr)
    {
        Pml<TGrid>::save(ostr);

        const int size = this->splitGrid->getNumPmlNodes();
        ostr.write((char*)&size, sizeof(size));

        ostr.write((char*)bCoeffX.data(), sizeof(FP) * size);
        ostr.write((char*)bCoeffY.data(), sizeof(FP) * size);
        ostr.write((char*)bCoeffZ.data(), sizeof(FP) * size);

        ostr.write((char*)eCoeffX.data(), sizeof(FP) * size);
        ostr.write((char*)eCoeffY.data(), sizeof(FP) * size);
        ostr.write((char*)eCoeffZ.data(), sizeof(FP) * size);
    }

    template<class TGrid>
    inline void PmlSpectral<TGrid>::load(std::istream& istr)
    {
        Pml<TGrid>::load(istr);

        int size = 0;
        istr.read((char*)&size, sizeof(size));

        bCoeffX.resize(size);
        bCoeffY.resize(size);
        bCoeffZ.resize(size);

        eCoeffX.resize(size);
        eCoeffY.resize(size);
        eCoeffZ.resize(size);

        istr.read((char*)bCoeffX.data(), sizeof(FP) * size);
        istr.read((char*)bCoeffY.data(), sizeof(FP) * size);
        istr.read((char*)bCoeffZ.data(), sizeof(FP) * size);

        istr.read((char*)eCoeffX.data(), sizeof(FP) * size);
        istr.read((char*)eCoeffY.data(), sizeof(FP) * size);
        istr.read((char*)eCoeffZ.data(), sizeof(FP) * size);
    }
}
