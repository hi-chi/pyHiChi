#pragma once
#include "FieldBoundaryCondition.h"

namespace pfc
{
    template <GridTypes gridTypes>
    class PeriodicalBoundaryConditionSpectral : public FieldBoundaryCondition<gridTypes> {
    public:

        // periodical boundaries for spectral solvers are default because of FFT
        // this class does nothing

        PeriodicalBoundaryConditionSpectral(CoordinateEnum axis,
            SpectralFieldSolver<gridTypes>* fieldSolver) :
            FieldBoundaryCondition<gridTypes>(axis, fieldSolver) {
        }

        void generateB() override {};
        void generateE() override {};

        FieldBoundaryCondition<gridTypes>* createInstance(
            FieldSolver<gridTypes>* fieldSolver = nullptr) override {
            SpectralFieldSolver<gridTypes>* solver = fieldSolver ?
                static_cast<SpectralFieldSolver<gridTypes>*>(fieldSolver) :
                static_cast<SpectralFieldSolver<gridTypes>*>(this->fieldSolver);
            return new PeriodicalBoundaryConditionSpectral(this->axis, solver);
        }
    };
}
