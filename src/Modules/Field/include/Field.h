#pragma once
#include "Vectors.h"
#include "AnalyticalFieldSolver.h"

#include <memory>


namespace pfc {

    // Field solver + computational grid
    template <class TFieldSolver>
    class Field {
    public:

        using TGrid = typename TFieldSolver::GridType;

        Field() {}
        Field(const Int3& numInternalCells,
            const FP3& minCoords, const FP3& steps, FP dt) {
            this->grid.reset(new TGrid(Int3(numInternalCells), minCoords, steps, numInternalCells));
            this->fieldSolver.reset(new TFieldSolver(grid.get(), dt));
            this->refresh();
        }

        void setTimeStep(FP dt) { this->fieldSolver->setTimeStep(dt); }
        FP getTimeStep() { return this->fieldSolver->dt; }

        void updateFields() { this->fieldSolver->updateFields(); }

        void setTime(FP t) { this->fieldSolver->globalTime = t; }
        FP getTime() { return this->fieldSolver->globalTime; }

        void refresh() {
            this->fieldSolver->globalTime = 0.0;
        }

        void advance(FP dt) {
            if (this->fieldSolver->dt != dt) {
                // all field solver submodules are refreshed
                this->fieldSolver->setTimeStep(dt);
            }
            this->fieldSolver->updateFields();
        }

        TGrid* getGrid() {
            return grid.get();
        }

        TFieldSolver* getFieldSolver() {
            return fieldSolver.get();
        }

        std::unique_ptr<TGrid> grid;
        std::unique_ptr<TFieldSolver> fieldSolver;

    };


    // Specialization for analytical field
    template<>
    class Field<AnalyticalFieldSolver> {
    public:

        Field() {}
        Field(FP dt) {
            this->fieldSolver.reset(new AnalyticalFieldSolver(dt));
            this->refresh();
        }

        void setTimeStep(FP dt) { this->fieldSolver->setTimeStep(dt); }
        FP getTimeStep() { return this->fieldSolver->dt; }

        void updateFields() { this->fieldSolver->updateFields(); }

        void setTime(FP t) { this->fieldSolver->field->globalTime = t; }
        FP getTime() { return this->fieldSolver->field->globalTime; }

        void refresh() {
            this->fieldSolver->field->globalTime = 0.0;
        }

        void advance(FP dt) {
            this->setTimeStep(dt);
            this->updateFields();
        }

        AnalyticalFieldSolver* getFieldSolver() {
            return this->fieldSolver.get();
        }

        AnalyticalField* getGrid() {
            return this->fieldSolver->field.get();
        }

        std::unique_ptr<AnalyticalFieldSolver> fieldSolver;

    };

}
