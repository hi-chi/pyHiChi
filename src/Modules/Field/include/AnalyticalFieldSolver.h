#pragma once
#include "AnalyticalField.h"


namespace pfc {

    // Contains general field solver methods
    class AnalyticalFieldSolver
    {
    public:

        using FunctionType = typename AnalyticalField::FunctionType;

        AnalyticalFieldSolver(FP dt) : field(new AnalyticalField()), dt(dt) {}

        AnalyticalFieldSolver(FP dt,
            FunctionType funcEx, FunctionType funcEy, FunctionType funcEz,
            FunctionType funcBx, FunctionType funcBy, FunctionType funcBz,
            FunctionType funcJx, FunctionType funcJy, FunctionType funcJz) : dt(dt)
        {
            field.reset(new AnalyticalField(funcEx, funcEy, funcEz,
                funcBx, funcBy, funcBz, funcJx, funcJy, funcJz));
        }

        AnalyticalFieldSolver(FP dt,
            FunctionType funcEx, FunctionType funcEy, FunctionType funcEz,
            FunctionType funcBx, FunctionType funcBy, FunctionType funcBz) : dt(dt)
        {
            field.reset(new AnalyticalField(funcEx, funcEy, funcEz, funcBx, funcBy, funcBz));
        }

        void updateFields() {
            this->field->globalTime += this->dt;
        }

        void setTimeStep(FP dt) {
            this->dt = dt;
        }

        std::unique_ptr<AnalyticalField> field;
        FP dt = 0.0;
    };
}
