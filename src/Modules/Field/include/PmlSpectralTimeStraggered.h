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
        PmlSpectralTimeStraggered(SpectralFieldSolver<gridTypes>* solver, Int3 sizePML) :
            PmlSpectral<gridTypes>((SpectralFieldSolver<gridTypes>*)solver, sizePML)
        {
            Int3 sizeComplex = FourierTransform::getSizeOfComplex(solver->grid->numCells);
            Int3 sizeStorage(sizeComplex.x, sizeComplex.y, 2 * sizeComplex.z);
            tmpFieldReal = ScalarField<FP>(sizeStorage);
            tmpFieldComplex = ScalarField<complexFP>(reinterpret_cast<complexFP*>(tmpFieldReal.getData()), sizeComplex);
        }

        virtual void updateBxSplit();
        virtual void updateBySplit();
        virtual void updateBzSplit();

        virtual void updateExSplit();
        virtual void updateEySplit();
        virtual void updateEzSplit();

        virtual void updateBSplit();
        virtual void updateESplit();

        virtual void doFirstStep();

        virtual void computeTmpField(MemberOfFP3 coordK, ScalarField<complexFP>& field, double dt) {};

        ScalarField<FP> tmpFieldReal;
        ScalarField<complexFP> tmpFieldComplex;  // has shared memory with tmpFieldReal

    protected:
        SpectralFieldSolver<gridTypes>* getFieldSolver() {
            return (SpectralFieldSolver<gridTypes>*)this->fieldSolver;
        }

        void updateFieldSplit(std::vector<FP>& field, double sign);
    };

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateFieldSplit(std::vector<FP>& field, double sign) {
#pragma omp parallel for
        for (int idx = 0; idx < this->numCells; ++idx)
        {
            int i = this->cellIndex[idx].x;
            int j = this->cellIndex[idx].y;
            int k = this->cellIndex[idx].z;

            field[idx] += sign * tmpFieldReal(i, j, k);

        }
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateBxSplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::y, fs->complexGrid->Ez, fs->grid->dt * 0.5);
        updateFieldSplit(this->byx, -1);
        computeTmpField(&FP3::z, fs->complexGrid->Ey, fs->grid->dt * 0.5);
        updateFieldSplit(this->bzx, +1);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateBySplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::z, fs->complexGrid->Ex, fs->grid->dt * 0.5);
        updateFieldSplit(this->bzy, -1);
        computeTmpField(&FP3::x, fs->complexGrid->Ez, fs->grid->dt * 0.5);
        updateFieldSplit(this->bxy, +1);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateBzSplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::x, fs->complexGrid->Ey, fs->grid->dt * 0.5);
        updateFieldSplit(this->bxz, -1);
        computeTmpField(&FP3::y, fs->complexGrid->Ex, fs->grid->dt * 0.5);
        updateFieldSplit(this->byz, +1);
    }


    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateExSplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::y, fs->complexGrid->Bz, fs->grid->dt);
        updateFieldSplit(this->eyx, +1);
        computeTmpField(&FP3::z, fs->complexGrid->By, fs->grid->dt);
        updateFieldSplit(this->ezx, -1);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateEySplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::z, fs->complexGrid->Bx, fs->grid->dt);
        updateFieldSplit(this->ezy, +1);
        computeTmpField(&FP3::x, fs->complexGrid->Bz, fs->grid->dt);
        updateFieldSplit(this->exy, -1);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateEzSplit()
    {
        SpectralFieldSolver<gridTypes>* fs = getFieldSolver();
        computeTmpField(&FP3::x, fs->complexGrid->By, fs->grid->dt);
        updateFieldSplit(this->exz, +1);
        computeTmpField(&FP3::y, fs->complexGrid->Bx, fs->grid->dt);
        updateFieldSplit(this->eyz, -1);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateBSplit()
    {
        updateBxSplit();
        updateBySplit();
        updateBzSplit();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::updateESplit()
    {
        updateExSplit();
        updateEySplit();
        updateEzSplit();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectralTimeStraggered<gridTypes>::doFirstStep()
    {
        updateBSplit();
        updateESplit();
        updateBSplit();
    }

}