#pragma once
#include "Vectors.h"
#include "AnalyticalField.h"


namespace pfc {

    // Field with a computational grid
    template <class TGrid, class TFieldSolver>
    class FieldEntity {
    public:

        FieldEntity() {}
        FieldEntity(const Int3& numInternalCells,
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

        void save(std::ostream& ostr) {
            ostr.write((char*)&(fieldSolver->dt), sizeof(fieldSolver->dt));
            grid->save(ostr);
            fieldSolver->save(ostr);
        }

        void load(std::istream& istr) {
            FP dt = 0.0;
            istr.read((char*)&dt, sizeof(dt));

            grid.reset(new TGrid());
            grid->load(istr);

            fieldSolver.reset(new TFieldSolver(grid.get(), dt));
            fieldSolver->load(istr);
        }

        std::unique_ptr<TGrid> grid;
        std::unique_ptr<TFieldSolver> fieldSolver;

    };

    class NoFieldSolver {};
    class NoGrid {};

    // Specialization for analytical field
    template<>
    class FieldEntity<NoGrid, NoFieldSolver> {
    public:

        FieldEntity() {}
        FieldEntity(FP dt) {
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

        void save(std::ostream& ostr) {
            field->save(ostr);
        }

        void load(std::istream& istr) {
            field->load(istr);
        }

        std::unique_ptr<AnalyticalField> field;

    };

}