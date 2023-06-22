#pragma once
#include "PmlSpectralTimeStaggered.h"

namespace pfc {

    template <class TGrid>
    class PmlPsatdTimeStaggered : public PmlSpectralTimeStaggered<TGrid, PmlPsatdTimeStaggered<TGrid>>
    {
    public:
        PmlPsatdTimeStaggered(TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt, Int3 sizePML) :
            PmlSpectralTimeStaggered<TGrid, PmlPsatdTimeStaggered<TGrid>>(grid, complexGrid, dt, sizePML) {}

        void computeTmpField(CoordinateEnum coordK,
            SpectralScalarField<FP, complexFP>& field, double dt);
    };

    template <class TGrid>
    inline void PmlPsatdTimeStaggered<TGrid>::computeTmpField(
        CoordinateEnum coordK, SpectralScalarField<FP, complexFP>& field, double dt)
    {
        // TODO: check border indices
        Int3 begin = Int3(0, 0, 0);
        Int3 end = this->complexGrid->numCells;

        OMP_FOR()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++) {
                    FP3 K = this->getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();
                    if (normK == 0) continue;
                    K = K / normK;

                    this->tmpFieldComplex(i, j, k) = 2.0 * sin(normK * constants::c * dt * 0.5) *
                        complexFP::i() * (complexFP)(K[(int)coordK]) * field(i, j, k);
                }

        this->fourierTransform.doFourierTransform(fourier_transform::Direction::CtoR);
    }
}
