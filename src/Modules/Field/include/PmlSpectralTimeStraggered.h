#pragma once
#include "Grid.h"
#include "FieldSolver.h"
#include "Pml.h"
#include "Constants.h"

namespace pfc {

    template<GridTypes gridTypes>
    class PmlSpectralTimeStraggered : public PmlSpectral<gridTypes>
    {
    public:
        PmlSpectralTimeStraggered(SpectralFieldSolver<gridTypes>* solver, Int3 sizePml) :
            PmlSpectral<gridTypes>(solver, sizePml),
            tmpFieldReal(solver->grid->sizeStorage),
            tmpFieldComplex(&tmpFieldReal, fourier_transform::getSizeOfComplexArray(solver->grid->numCells))
        {
            fourierTransform.initialize(&tmpFieldReal, &tmpFieldComplex, solver->grid->numCells);
        }

        void updateBxSplit() override;
        void updateBySplit() override;
        void updateBzSplit() override;

        void updateExSplit() override;
        void updateEySplit() override;
        void updateEzSplit() override;

        virtual void computeTmpField(MemberOfFP3 coordK,
            SpectralScalarField<FP, complexFP>& field, double dt, double sign) {};

        ScalarField<FP> tmpFieldReal;
        SpectralScalarField<FP, complexFP> tmpFieldComplex;
        FourierTransformField fourierTransform;

    protected:
        SpectralFieldSolver<gridTypes>* getFieldSolver() {
            return static_cast<SpectralFieldSolver<gridTypes>*>(this->fieldSolver);
        }

        void updateFieldSplit(std::vector<FP>& field,
            const std::vector<Int3>& index);
    };

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateFieldSplit(std::vector<FP>& field,
        const std::vector<Int3>& index)
    {
        const int size = index.size();

        OMP_FOR()
            for (int idx = 0; idx < size; ++idx)
            {
                int i = index[idx].x, j = index[idx].y, k = index[idx].z;

                field[idx] += tmpFieldReal(i, j, k);
            }
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateBxSplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::y, fs->complexGrid->Ez, fs->dt * 0.5, -1);
        updateFieldSplit(this->byx, this->bIndex);
        computeTmpField(&FP3::z, fs->complexGrid->Ey, fs->dt * 0.5, +1);
        updateFieldSplit(this->bzx, this->bIndex);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateBySplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::z, fs->complexGrid->Ex, fs->dt * 0.5, -1);
        updateFieldSplit(this->bzy, this->bIndex);
        computeTmpField(&FP3::x, fs->complexGrid->Ez, fs->dt * 0.5, +1);
        updateFieldSplit(this->bxy, this->bIndex);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateBzSplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::x, fs->complexGrid->Ey, fs->dt * 0.5, -1);
        updateFieldSplit(this->bxz, this->bIndex);
        computeTmpField(&FP3::y, fs->complexGrid->Ex, fs->dt * 0.5, +1);
        updateFieldSplit(this->byz, this->bIndex);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateExSplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::y, fs->complexGrid->Bz, fs->dt, +1);
        updateFieldSplit(this->eyx, this->eIndex);
        computeTmpField(&FP3::z, fs->complexGrid->By, fs->dt, -1);
        updateFieldSplit(this->ezx, this->eIndex);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateEySplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::z, fs->complexGrid->Bx, fs->dt, +1);
        updateFieldSplit(this->ezy, this->eIndex);
        computeTmpField(&FP3::x, fs->complexGrid->Bz, fs->dt, -1);
        updateFieldSplit(this->exy, this->eIndex);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateEzSplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::x, fs->complexGrid->By, fs->dt, +1);
        updateFieldSplit(this->exz, this->eIndex);
        computeTmpField(&FP3::y, fs->complexGrid->Bx, fs->dt, -1);
        updateFieldSplit(this->eyz, this->eIndex);
    }

}
