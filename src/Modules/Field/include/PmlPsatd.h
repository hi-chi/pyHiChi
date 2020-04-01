#pragma once
#include "Grid.h"
#include "FieldSolver.h"
#include "PmlSpectralTimeStraggered.h"
#include "Constants.h"

namespace pfc {

    class PmlPsatd : public PmlSpectralTimeStraggered<GridTypes::PSATDGridType>
    {
    public:
        PmlPsatd(SpectralFieldSolver<GridTypes::PSATDGridType>* solver, Int3 sizePML) :
            PmlSpectralTimeStraggered((SpectralFieldSolver<GridTypes::PSATDGridType>*)solver, sizePML) {}
        
        virtual void computeTmpField(MemberOfFP3 coordK, ScalarField<complexFP>& field, double dt);
    };

    inline void PmlPsatd::computeTmpField(MemberOfFP3 coordK, ScalarField<complexFP>& field, double dt)
    {
        SpectralFieldSolver<GridTypes::PSATDGridType>* fs = getFieldSolver();
        Int3 begin = fs->updateComplexBAreaBegin;
        Int3 end = fs->updateComplexBAreaEnd;
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++) {
                    FP3 K = fs->getWaveVector(Int3(i, j, k));
                    FP normK = K.norm();
                    if (normK == 0) continue;
                    K = K / normK;

                    tmpFieldComplex(i, j, k) = 2 * sin(normK*constants::c*dt*0.5) * complexFP::i() *
                        (complexFP)(K.*coordK) * field(i, j, k);
                }
        FourierTransform::doFourierTransform(tmpFieldReal, tmpFieldComplex, FourierTransformDirection::CtoR);
    }

}