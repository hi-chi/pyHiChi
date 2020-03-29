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
        PmlPsatdBase(SpectralFieldSolver<TPSATDGridType>* solver, Int3 sizePML) :
            PmlSpectralTimeStraggered<TPSATDGridType>((SpectralFieldSolver<TPSATDGridType>*)solver, sizePML) {}

        void computeTmpField(MemberOfFP3 coordK, ScalarField<complexFP>& field, double dt);
    };

    template <GridTypes TPSATDGridType>
    inline void PmlPsatdBase<TPSATDGridType>::computeTmpField(MemberOfFP3 coordK, ScalarField<complexFP>& field, double dt)
    {
        SpectralFieldSolver<TPSATDGridType>* fs = PmlSpectralTimeStraggered<TPSATDGridType>::getFieldSolver();
        Int3 begin = fs->updateComplexBAreaBegin;
        Int3 end = fs->updateComplexBAreaEnd;
#pragma omp parallel for
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
        FourierTransform::doFourierTransform(PmlSpectralTimeStraggered<TPSATDGridType>::tmpFieldReal,
            PmlSpectralTimeStraggered<TPSATDGridType>::tmpFieldComplex, FourierTransformDirection::CtoR);
    }


    class PmlPsatdTimeStraggered : public PmlPsatdBase<GridTypes::PSATDTimeStraggeredGridType>
    {
    public:
        PmlPsatdTimeStraggered(SpectralFieldSolver<GridTypes::PSATDTimeStraggeredGridType>* solver, Int3 sizePML) :
            PmlPsatdBase<GridTypes::PSATDTimeStraggeredGridType>((SpectralFieldSolver<GridTypes::PSATDTimeStraggeredGridType>*)solver, sizePML) {}
        
        virtual void computeTmpField(MemberOfFP3 coordK, ScalarField<complexFP>& field, double dt) {
            PmlPsatdBase<GridTypes::PSATDTimeStraggeredGridType>::computeTmpField(coordK, field, dt);
        }
    };


    class PmlPsatd : public PmlPsatdBase<GridTypes::PSATDGridType>
    {
    public:
        PmlPsatd(SpectralFieldSolver<GridTypes::PSATDGridType>* solver, Int3 sizePML) :
            PmlPsatdBase<GridTypes::PSATDGridType>((SpectralFieldSolver<GridTypes::PSATDGridType>*)solver, sizePML) {}

        virtual void computeTmpField(MemberOfFP3 coordK, ScalarField<complexFP>& field, double dt) {
            PmlPsatdBase<GridTypes::PSATDGridType>::computeTmpField(coordK, field, dt);
        }
    };

}