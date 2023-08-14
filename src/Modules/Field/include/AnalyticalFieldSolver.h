#pragma once
#include "AnalyticalField.h"


namespace pfc {

    // Contains general field solver methods
    class AnalyticalFieldSolver
    {
    public:

        using GridType = AnalyticalField;  // AnalyticalField plays a grid role

        AnalyticalFieldSolver(AnalyticalField* field, FP dt) : field(field), dt(dt) {}

        // constructor for loading
        explicit AnalyticalFieldSolver(AnalyticalField* field) : field(field) {}

        void updateFields() {
            this->field->globalTime += this->dt;
        }

        void setTimeStep(FP dt) {
            this->dt = dt;
        }
        FP getTimeStep() const {
            return this->dt;
        }

        void setTime(FP t) {
            this->field->globalTime = t;
        }
        FP getTime() const {
            return this->field->globalTime;
        }

        void save(std::ostream& ostr) {
            ostr.write((char*)&dt, sizeof(dt));
        }
        void load(std::istream& istr) {
            istr.read((char*)&dt, sizeof(dt));
        }

        AnalyticalField* field = nullptr;
        FP dt = 0.0;
    };
}
