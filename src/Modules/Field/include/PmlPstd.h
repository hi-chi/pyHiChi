#pragma once
#include "PmlSpectralTimeStaggered.h"

namespace pfc {

    class PmlPstd : public PmlSpectralTimeStaggered<PSTDGrid, PmlPstd>
    {
    public:
        PmlPstd(PSTDGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt, Int3 sizePML) :
            PmlSpectralTimeStaggered(grid, complexGrid, dt, sizePML) {}

        void computeTmpField(CoordinateEnum coordK,
            SpectralScalarField<FP, complexFP>& field, double dt);
    };

    inline void PmlPstd::computeTmpField(CoordinateEnum coordK,
        SpectralScalarField<FP, complexFP>& field, double dt)
    {
        // TODO: check border indices
        Int3 begin = Int3(0, 0, 0);
        Int3 end = this->complexGrid->numCells;

        OMP_FOR()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++) {
                    this->tmpFieldComplex(i, j, k) = constants::c * dt * complexFP::i() *
                        (complexFP)(this->getWaveVector(Int3(i, j, k))[(int)coordK]) * field(i, j, k);
                }

        this->fourierTransform.doFourierTransform(fourier_transform::Direction::CtoR);
    }

}
