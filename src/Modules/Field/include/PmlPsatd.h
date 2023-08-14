#pragma once
#include "PmlSpectralTimeStaggered.h"

namespace pfc {

    template <class TGrid>
    class PmlPsatdTimeStaggered : public PmlSpectralTimeStaggered<TGrid, PmlPsatdTimeStaggered<TGrid>>
    {
    public:
        PmlPsatdTimeStaggered(TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt,
            Int3 domainIndexBegin, Int3 domainIndexEnd, Int3 complexDomainIndexBegin, Int3 complexDomainIndexEnd,
            Int3 sizePML, FP nPmlParam = (FP)4.0, FP r0PmlParam = (FP)1e-8) :
            PmlSpectralTimeStaggered<TGrid, PmlPsatdTimeStaggered<TGrid>>(grid, complexGrid, dt,
                domainIndexBegin, domainIndexEnd, complexDomainIndexBegin, complexDomainIndexEnd,
                sizePML, nPmlParam, r0PmlParam)
        {}

        // constructor for loading
        explicit PmlPsatdTimeStaggered(TGrid* grid, SpectralGrid<FP, complexFP>* complexGrid, FP dt,
            Int3 domainIndexBegin, Int3 domainIndexEnd, Int3 complexDomainIndexBegin, Int3 complexDomainIndexEnd) :
            PmlSpectralTimeStaggered<TGrid, PmlPsatdTimeStaggered<TGrid>>(grid, complexGrid, dt,
                domainIndexBegin, domainIndexEnd, complexDomainIndexBegin, complexDomainIndexEnd)
        {}

        void computeTmpField(CoordinateEnum coordK,
            SpectralScalarField<FP, complexFP>& field, double dt);
    };

    template <class TGrid>
    inline void PmlPsatdTimeStaggered<TGrid>::computeTmpField(
        CoordinateEnum coordK, SpectralScalarField<FP, complexFP>& field, double dt)
    {
        const Int3 begin = this->complexDomainIndexBegin;
        const Int3 end = this->complexDomainIndexEnd;

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
