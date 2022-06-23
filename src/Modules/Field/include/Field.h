#pragma once
#include "Vectors.h"
#include "AnalyticalField.h"


namespace pfc {

    // Field solver + computational grid
    template <class TGrid, class TFieldSolver>
    class Field {
    public:

        Field() {}
        Field(const Int3& numInternalCells,
            const FP3& minCoords, const FP3& steps, FP dt) {
            grid.reset(new TGrid(Int3(numInternalCells), minCoords, steps, numInternalCells));
            fieldSolver.reset(new TFieldSolver(grid.get(), dt));
        }

        void refresh() {
            this->fieldSolver->globalTime = 0.0;
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


    class NoFieldSolver {};
    class NoGrid {};

    // Specialization for analytical field
    template<>
    class Field<NoGrid, NoFieldSolver> {
    public:

        Field() {}
        Field(FP dt) {
            field.reset(new AnalyticalField(dt));
        }

        void refresh() {
            this->field->globalTime = 0.0;
        }

        AnalyticalField* getAnalyticalField() {
            return field.get();
        }

        AnalyticalField* getGrid() {
            return field.get();
        }

        AnalyticalField* getFieldSolver() {
            return field.get();
        }

        std::unique_ptr<AnalyticalField> field;

    };

}
