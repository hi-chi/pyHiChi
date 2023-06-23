#pragma once
#include "Pml.h"

namespace pfc {

    template<class TGrid, class TDerived>
    class PmlSpectralTimeStaggered : public PmlSpectral<TGrid>
    {
    public:
        PmlSpectralTimeStaggered(TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt, Int3 sizePML,
            Int3 domainIndexBegin, Int3 domainIndexEnd, Int3 complexDomainIndexBegin, Int3 complexDomainIndexEnd) :
            PmlSpectral<TGrid>(grid, complexGrid, dt, sizePML, domainIndexBegin, domainIndexEnd,
                complexDomainIndexBegin, complexDomainIndexEnd),
            tmpFieldReal(this->grid->sizeStorage),
            tmpFieldComplex(&tmpFieldReal, fourier_transform::getSizeOfComplexArray(this->grid->numCells))
        {
            fourierTransform.initialize(&tmpFieldReal, &tmpFieldComplex, this->grid->numCells);
        }

        /* implement the next methods in derived classes
        void computeTmpField(CoordinateEnum coordK, SpectralScalarField<FP, complexFP>& field, double dt);
        */

        void updateBSplit();
        void updateESplit();

        FP3 getWaveVector(const Int3& ind);

        ScalarField<FP> tmpFieldReal;
        SpectralScalarField<FP, complexFP> tmpFieldComplex;
        FourierTransformField fourierTransform;
    };

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateBSplit()
    {
        TDerived* derived = static_cast<TDerived*>(this);

        const int size = this->splitGrid->getNumPmlNodes();

        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Ez, this->dt * 0.5);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->bxy[idx] -= this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->Ey, this->dt * 0.5);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->bxz[idx] += this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->Ex, this->dt * 0.5);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->byz[idx] -= this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->Ez, this->dt * 0.5);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->byx[idx] += this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->Ey, this->dt * 0.5);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->bzx[idx] -= this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Ex, this->dt * 0.5);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->bzy[idx] += this->tmpFieldReal(this->splitGrid->getIndex3d(idx));
    }

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateESplit()
    {
        TDerived* derived = static_cast<TDerived*>(this);

        const int size = this->splitGrid->getNumPmlNodes();

        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Bz, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->exy[idx] += this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->By, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->exz[idx] -= this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->Bx, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->eyz[idx] += this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->Bz, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->eyx[idx] -= this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->By, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->ezx[idx] += this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Bx, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->ezy[idx] -= this->tmpFieldReal(this->splitGrid->getIndex3d(idx));
    }

    template<class TGrid, class TDerived>
    inline FP3 PmlSpectralTimeStaggered<TGrid, TDerived>::getWaveVector(const Int3& ind)
    {
        // TODO: check number of cells
        FP kx = (2 * constants::pi * ((ind.x <= this->grid->numCells.x / 2) ? ind.x : ind.x - this->grid->numCells.x)) /
            (this->grid->steps.x * this->grid->numCells.x);
        FP ky = (2 * constants::pi * ((ind.y <= this->grid->numCells.y / 2) ? ind.y : ind.y - this->grid->numCells.y)) /
            (this->grid->steps.y * this->grid->numCells.y);
        FP kz = (2 * constants::pi * ((ind.z <= this->grid->numCells.z / 2) ? ind.z : ind.z - this->grid->numCells.z)) /
            (this->grid->steps.z * this->grid->numCells.z);
        return FP3(kx, ky, kz);
    }

}
