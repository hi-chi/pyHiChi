#pragma once
#include "Grid.h"
#include "FieldSolver.h"
#include "PmlSpectralTimeStraggered.h"
#include "Constants.h"

namespace pfc {

    template <GridTypes TPSATDGridType>
    class PmlPsatdBase : public PmlSpectralTimeStraggered<TPSATDGridType>
    {
    public:
        PmlPsatdBase(SpectralFieldSolver<TPSATDGridType>* solver, Int3 sizePml) :
            PmlSpectralTimeStraggered<TPSATDGridType>(solver, sizePml) {}

        void computeTmpField(MemberOfFP3 coordK,
            SpectralScalarField<FP, complexFP>& field, double dt, double sign);
    };

    template <GridTypes TPSATDGridType>
    inline void PmlPsatdBase<TPSATDGridType>::computeTmpField(MemberOfFP3 coordK,
        SpectralScalarField<FP, complexFP>& field, double dt, double sign)
    {
        SpectralFieldSolver<TPSATDGridType>* fs = PmlSpectralTimeStraggered<TPSATDGridType>::getFieldSolver();
        const Int3 begin = fs->updateComplexBAreaBegin;
        const Int3 end = fs->updateComplexBAreaEnd;

        OMP_FOR_COLLAPSE()
            for (int i = begin.x; i < end.x; i++)
                for (int j = begin.y; j < end.y; j++)
                    for (int k = begin.z; k < end.z; k++) {
                        FP3 K = fs->getWaveVector(Int3(i, j, k));
                        FP normK = K.norm();
                        if (normK == 0) continue;
                        K = K / normK;

                        PmlSpectralTimeStraggered<TPSATDGridType>::tmpFieldComplex(i, j, k) =
                            sign * 2.0 * sin(normK * constants::c * dt * 0.5) * complexFP::i() *
                            (complexFP)(K.*coordK) * field(i, j, k);
                    }

        PmlSpectralTimeStraggered<TPSATDGridType>::fourierTransform.doFourierTransform(fourier_transform::Direction::CtoR);
    }


    class PmlPsatdTimeStraggered : public PmlPsatdBase<GridTypes::PSATDTimeStraggeredGridType>
    {
    public:
        PmlPsatdTimeStraggered(SpectralFieldSolver<GridTypes::PSATDTimeStraggeredGridType>* solver, Int3 sizePml) :
            PmlPsatdBase<GridTypes::PSATDTimeStraggeredGridType>(solver, sizePml) {}

        virtual void computeTmpField(MemberOfFP3 coordK,
            SpectralScalarField<FP, complexFP>& field, double dt, double sign) {
            PmlPsatdBase<GridTypes::PSATDTimeStraggeredGridType>::computeTmpField(coordK, field, dt, sign);
        }
    };


    class PmlPsatd : public PmlPsatdBase<GridTypes::PSATDGridType>
    {
    public:
        PmlPsatd(SpectralFieldSolver<GridTypes::PSATDGridType>* solver, Int3 sizePml) :
            PmlPsatdBase<GridTypes::PSATDGridType>(solver, sizePml) {}

        virtual void computeTmpField(MemberOfFP3 coordK,
            SpectralScalarField<FP, complexFP>& field, double dt, double sign) {
            PmlPsatdBase<GridTypes::PSATDGridType>::computeTmpField(coordK, field, dt, sign);
        }
    };

}
