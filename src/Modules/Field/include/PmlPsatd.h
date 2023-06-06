#pragma once
#include "Grid.h"
#include "FieldSolver.h"
#include "PmlSpectralTimeStraggered.h"
#include "Constants.h"
#include "macros.h"

namespace pfc {

    template <GridTypes TPSATDGridType>
    class PmlPsatdTimeStraggered : public PmlSpectralTimeStraggered<TPSATDGridType>
    {
    public:
        PmlPsatdTimeStraggered(SpectralFieldSolver<TPSATDGridType>* solver, Int3 sizePML) :
            PmlSpectralTimeStraggered<TPSATDGridType>(solver, sizePML) {}

        void computeTmpField(MemberOfFP3 coordK,
            SpectralScalarField<FP, complexFP>& field, double dt);
    };

    template <GridTypes TPSATDGridType>
    inline void PmlPsatdTimeStraggered<TPSATDGridType>::computeTmpField(MemberOfFP3 coordK,
        SpectralScalarField<FP, complexFP>& field, double dt)
    {
        SpectralFieldSolver<TPSATDGridType>* fs = PmlSpectralTimeStraggered<TPSATDGridType>::getFieldSolver();
        Int3 begin = fs->updateComplexBAreaBegin;
        Int3 end = fs->updateComplexBAreaEnd;

        OMP_FOR_COLLAPSE()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++) {
                    FP3 K = fs->getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();
                    if (normK == 0) continue;
                    K = K / normK;

                    PmlSpectralTimeStraggered<TPSATDGridType>::tmpFieldComplex(i, j, k) =
                        2.0 * sin(normK*constants::c*dt*0.5) * complexFP::i() *
                        (complexFP)(K.*coordK) * field(i, j, k);
                }
        PmlSpectralTimeStraggered<TPSATDGridType>::fourierTransform.doFourierTransform(fourier_transform::Direction::CtoR);
    }
}