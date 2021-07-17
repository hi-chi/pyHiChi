#pragma once
#include "Grid.h"
#include "FieldSolver.h"
#include "PmlSpectralTimeStraggered.h"
#include "Constants.h"

namespace pfc {

    class PmlPstd : public PmlSpectralTimeStraggered<GridTypes::PSTDGridType>
    {
    public:
        PmlPstd(SpectralFieldSolver<GridTypes::PSTDGridType>* solver, Int3 sizePML) :
            PmlSpectralTimeStraggered((SpectralFieldSolver<GridTypes::PSTDGridType>*)solver, sizePML) {}

        virtual void computeTmpField(MemberOfFP3 coordK,
            SpectralScalarField<FP, complexFP>& field, double dt);
    };

    inline void PmlPstd::computeTmpField(MemberOfFP3 coordK,
        SpectralScalarField<FP, complexFP>& field, double dt)
    {
        SpectralFieldSolver<GridTypes::PSTDGridType>* fs = getFieldSolver();
        Int3 begin = fs->updateComplexBAreaBegin;
        Int3 end = fs->updateComplexBAreaEnd;
        OMP_FOR()
        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                    tmpFieldComplex(i, j, k) = constants::c * dt * complexFP::i() *
                    (complexFP)(fs->getWaveVector(Int3(i, j, k)).*coordK) * field(i, j, k);
        fourierTransform.doFourierTransform(fourier_transform::Direction::CtoR);
    }

}