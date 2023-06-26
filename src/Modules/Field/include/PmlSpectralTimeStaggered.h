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
            tmpFieldComplex(&tmpFieldReal, fourier_transform::getSizeOfComplexArray(domainIndexEnd - domainIndexBegin))
        {
            fourierTransform.initialize(&tmpFieldReal, &tmpFieldComplex, domainIndexEnd - domainIndexBegin);
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

        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Ez, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->bxy[idx] -= this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->Ey, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->bxz[idx] += this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->Ex, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->byz[idx] -= this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->Ez, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->byx[idx] += this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->Ey, this->dt);
        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
            this->splitGrid->bzx[idx] -= this->tmpFieldReal(this->splitGrid->getIndex3d(idx));

        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Ex, this->dt);
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
        const Int3 domainSize = this->domainIndexEnd - this->domainIndexBegin;
        FP kx = ((FP)2.0 * constants::pi * ((ind.x <= domainSize.x / 2) ? ind.x : ind.x - domainSize.x)) /
            (this->grid->steps.x * domainSize.x);
        FP ky = ((FP)2.0 * constants::pi * ((ind.y <= domainSize.y / 2) ? ind.y : ind.y - domainSize.y)) /
            (this->grid->steps.y * domainSize.y);
        FP kz = ((FP)2.0 * constants::pi * ((ind.z <= domainSize.z / 2) ? ind.z : ind.z - domainSize.z)) /
            (this->grid->steps.z * domainSize.z);
        return FP3(kx, ky, kz);
    }

}
