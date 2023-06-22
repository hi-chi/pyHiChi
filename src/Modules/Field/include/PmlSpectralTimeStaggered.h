#pragma once
#include "Pml.h"

namespace pfc {

    template<class TGrid, class TDerived>
    class PmlSpectralTimeStaggered : public PmlSpectral<TGrid>
    {
    public:
        PmlSpectralTimeStaggered(TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt, Int3 sizePML) :
            PmlSpectral<TGrid>(grid, complexGrid, dt, sizePML),
            tmpFieldReal(this->grid->sizeStorage),
            tmpFieldComplex(&tmpFieldReal, fourier_transform::getSizeOfComplexArray(this->grid->numCells))
        {
            fourierTransform.initialize(&tmpFieldReal, &tmpFieldComplex, this->grid->numCells);
        }

        /* implement the next methods in derived classes
        void computeTmpField(CoordinateEnum coordK, SpectralScalarField<FP, complexFP>& field, double dt);
        */
        
        void updateBxSplit();
        void updateBySplit();
        void updateBzSplit();

        void updateExSplit();
        void updateEySplit();
        void updateEzSplit();

        void updateBSplit();
        void updateESplit();

        FP3 getWaveVector(const Int3& ind);

        ScalarField<FP> tmpFieldReal;
        SpectralScalarField<FP, complexFP> tmpFieldComplex;
        FourierTransformField fourierTransform;

    protected:
        void updateFieldSplit(std::vector<FP>& field,
            const std::vector<Int3>& index, FP sign);
    };

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateFieldSplit(
        std::vector<FP>& field, const std::vector<Int3>& index, FP sign)
    {
        const int size = index.size();

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            int i = index[idx].x, j = index[idx].y, k = index[idx].z;

            field[idx] += sign * tmpFieldReal(i, j, k);
        }
    }

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateBxSplit()
    {
        TDerived* derived = static_cast<TDerived*>(this);
        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Ez, this->dt * 0.5);
        updateFieldSplit(this->byx, this->bIndex, -1);
        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->Ey, this->dt * 0.5);
        updateFieldSplit(this->bzx, this->bIndex, +1);
    }

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateBySplit()
    {
        TDerived* derived = static_cast<TDerived*>(this);
        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->Ex, this->dt * 0.5);
        updateFieldSplit(this->bzy, this->bIndex, -1);
        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->Ez, this->dt * 0.5);
        updateFieldSplit(this->bxy, this->bIndex, +1);
    }

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateBzSplit()
    {
        TDerived* derived = static_cast<TDerived*>(this);
        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->Ey, this->dt * 0.5);
        updateFieldSplit(this->bxz, this->bIndex, -1);
        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Ex, this->dt * 0.5);
        updateFieldSplit(this->byz, this->bIndex, +1);
    }

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateExSplit()
    {
        TDerived* derived = static_cast<TDerived*>(this);
        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Bz, this->dt);
        updateFieldSplit(this->eyx, this->eIndex, +1);
        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->By, this->dt);
        updateFieldSplit(this->ezx, this->eIndex, -1);
    }

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateEySplit()
    {
        TDerived* derived = static_cast<TDerived*>(this);
        derived->computeTmpField(CoordinateEnum::z, this->complexGrid->Bx, this->dt);
        updateFieldSplit(this->ezy, this->eIndex, +1);
        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->Bz, this->dt);
        updateFieldSplit(this->exy, this->eIndex, -1);
    }

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateEzSplit()
    {
        TDerived* derived = static_cast<TDerived*>(this);
        derived->computeTmpField(CoordinateEnum::x, this->complexGrid->By, this->dt);
        updateFieldSplit(this->exz, this->eIndex, +1);
        derived->computeTmpField(CoordinateEnum::y, this->complexGrid->Bx, this->dt);
        updateFieldSplit(this->eyz, this->eIndex, -1);
    }

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateBSplit()
    {
        updateBxSplit();
        updateBySplit();
        updateBzSplit();
    }

    template<class TGrid, class TDerived>
    inline void PmlSpectralTimeStaggered<TGrid, TDerived>::updateESplit()
    {
        updateExSplit();
        updateEySplit();
        updateEzSplit();
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
