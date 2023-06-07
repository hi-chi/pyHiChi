#pragma once
#include <memory>

#include "Grid.h"
#include "AnalyticalFieldSolver.h"
#include "FieldValue.h"
#include "Mapping.h"
#include "Fdtd.h"
#include "Psatd.h"
#include "PsatdTimeStaggered.h"
#include "Pstd.h"
#include "Mapping.h"
#include "Field.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace pfc
{
    using FunctionPointer = int64_t;

    // Some field solver properties
    template <class TFieldSolver>
    struct isFieldSolverSpatialStaggered {
        static const bool value = TFieldSolver::GridType::ifFieldsSpatialStraggered;
    };

    template <class TFieldSolver>
    struct isFieldSolverSupportPeriodicalBoundaryCondition {
        static const bool value = true;  // all solvers support periodical boundary conditions
    };

    template <class TFieldSolver>
    struct isFieldSolverSupportReflectBoundaryCondition {
        static const bool value = std::is_same<TFieldSolver, FDTD>::value;
    };

    template <class TFieldSolver>
    struct isFieldSolverSupportFieldGenerator {
        static const bool value = std::is_same<TFieldSolver, FDTD>::value;
    };

    template <class TFieldSolver>
    struct isFieldSolverSupportPml {
        static const bool value = true;  // all solvers support pml
    };

    template <class TFieldSolver>
    struct isFieldSolverSupportPoissonEquation {
        static const bool value = std::is_same<TFieldSolver, PSATD>::value ||
            std::is_same<TFieldSolver, PSATDTimeStraggered>::value ||
            std::is_same<TFieldSolver, PSATDPoisson>::value ||
            std::is_same<TFieldSolver, PSATDTimeStraggeredPoisson>::value;
    };

    template <class TFieldSolver>
    struct isFieldAnalytical {
        static const bool value = std::is_same<TFieldSolver, AnalyticalFieldSolver>::value;
    };


    // Interface depending on spatial template (shifted or collocated)
    template <class TFieldSolver, class TPyField, bool ifSpatialStraggered>
    class pySpatialStraggeredFieldInterface {};

    // Interface for spatial straggered grids
    template <class TFieldSolver, class TPyField>
    class pySpatialStraggeredFieldInterface<TFieldSolver, TPyField, true>
    {
    public:

        using TGrid = typename TFieldSolver::GridType;

        template <class FieldConfigurationType>
        void setFieldConfiguration(const FieldConfigurationType* fieldConf) {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            const int chunkSize = 32;
            const int nChunks = grid->numCells.z / chunkSize;
            const int chunkRem = grid->numCells.z % chunkSize;
            const int nx = grid->numCells.x, ny = grid->numCells.y;

            OMP_FOR_COLLAPSE()
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 cEx[chunkSize], cEy[chunkSize], cEz[chunkSize];
                        FP3 cBx[chunkSize], cBy[chunkSize], cBz[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
                        OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            cEx[k] = derived->convertCoords(grid->ExPosition(i, j, chunk * chunkSize));
                            cEy[k] = derived->convertCoords(grid->EyPosition(i, j, chunk * chunkSize));
                            cEz[k] = derived->convertCoords(grid->EzPosition(i, j, chunk * chunkSize));

                            cBx[k] = derived->convertCoords(grid->BxPosition(i, j, chunk * chunkSize));
                            cBy[k] = derived->convertCoords(grid->ByPosition(i, j, chunk * chunkSize));
                            cBz[k] = derived->convertCoords(grid->BzPosition(i, j, chunk * chunkSize));
                        }
                        OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            grid->Ex(i, j, k) = fieldConf->getE(cEx[k].x, cEx[k].y, cEx[k].z).x;
                            grid->Ey(i, j, k) = fieldConf->getE(cEy[k].x, cEy[k].y, cEy[k].z).y;
                            grid->Ez(i, j, k) = fieldConf->getE(cEz[k].x, cEz[k].y, cEz[k].z).z;

                            grid->Bx(i, j, k) = fieldConf->getB(cBx[k].x, cBx[k].y, cBx[k].z).x;
                            grid->By(i, j, k) = fieldConf->getB(cBy[k].x, cBy[k].y, cBy[k].z).y;
                            grid->Bz(i, j, k) = fieldConf->getB(cBz[k].x, cBz[k].y, cBz[k].z).z;
                        }
                    }
        }

    };

    // Interface for collocated grids
    template <class TFieldSolver, class TPyField>
    class pySpatialStraggeredFieldInterface<TFieldSolver, TPyField, false>
    {
    public:

        using TGrid = typename TFieldSolver::GridType;

        template <class FieldConfigurationType>
        void setFieldConfiguration(const FieldConfigurationType* fieldConf) {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            const int chunkSize = 32;
            const int nChunks = grid->numCells.z / chunkSize;
            const int chunkRem = grid->numCells.z % chunkSize;
            const int nx = grid->numCells.x, ny = grid->numCells.y;

            OMP_FOR_COLLAPSE()
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 coords[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
                        FP3 startPosition = grid->ExPosition(i, j, chunk * chunkSize);
                        OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            FP3 position(startPosition.x, startPosition.y, startPosition.z + k * grid->steps.z);
                            coords[k] = derived->convertCoords(position);
                        }
                        OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            FP3 E, B;
                            fieldConf->getEB(coords[k].x, coords[k].y, coords[k].z, &E, &B);

                            grid->Ex(i, j, k + chunk * chunkSize) = E.x;
                            grid->Ey(i, j, k + chunk * chunkSize) = E.y;
                            grid->Ez(i, j, k + chunk * chunkSize) = E.z;

                            grid->Bx(i, j, k + chunk * chunkSize) = B.x;
                            grid->By(i, j, k + chunk * chunkSize) = B.y;
                            grid->Bz(i, j, k + chunk * chunkSize) = B.z;
                        }
                    }
        }

        void pySetEMField(py::function fValueField)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            for (int i = 0; i < grid->numCells.x; i++)
                for (int j = 0; j < grid->numCells.y; j++)
                    for (int k = 0; k < grid->numCells.z; k++)
                    {
                        FP3 coords = derived->convertCoords(grid->ExPosition(i, j, k));
                        ValueField field;
                        fValueField("x"_a = coords.x, "y"_a = coords.y, "z"_a = coords.z,
                            "field_value"_a = std::reference_wrapper<ValueField>(field));

                        grid->Ex(i, j, k) = field.E.x;
                        grid->Ey(i, j, k) = field.E.y;
                        grid->Ez(i, j, k) = field.E.z;

                        grid->Bx(i, j, k) = field.B.x;
                        grid->By(i, j, k) = field.B.y;
                        grid->Bz(i, j, k) = field.B.z;
                    }
        }

        void setEMField(FunctionPointer _fValueField)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            void(*fValueField)(FP, FP, FP, FP*) = (void(*)(FP, FP, FP, FP*))_fValueField;
            const int chunkSize = 32;
            const int nChunks = grid->numCells.z / chunkSize;
            const int chunkRem = grid->numCells.z % chunkSize;
            const int nx = grid->numCells.x, ny = grid->numCells.y;

            OMP_FOR_COLLAPSE()
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 coords[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
                        FP3 startPosition = grid->ExPosition(i, j, chunk * chunkSize);
                        OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            FP3 position(startPosition.x, startPosition.y,
                                startPosition.z + k * grid->steps.z);
                            coords[k] = derived->convertCoords(position);
                        }
                        OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            ValueField field(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                            fValueField(coords[k].x, coords[k].y, coords[k].z, &(field.E.x));

                            grid->Ex(i, j, k + chunk * chunkSize) = field.E.x;
                            grid->Ey(i, j, k + chunk * chunkSize) = field.E.y;
                            grid->Ez(i, j, k + chunk * chunkSize) = field.E.z;

                            grid->Bx(i, j, k + chunk * chunkSize) = field.B.x;
                            grid->By(i, j, k + chunk * chunkSize) = field.B.y;
                            grid->Bz(i, j, k + chunk * chunkSize) = field.B.z;
                        }
                    }
        }

        void pyApplyFunction(py::function func)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            for (int i = 0; i < grid->numCells.x; i++)
                for (int j = 0; j < grid->numCells.y; j++)
                    for (int k = 0; k < grid->numCells.z; k++)
                    {
                        FP3 coords = derived->convertCoords(grid->ExPosition(i, j, k));
                        ValueField field(grid->Ex(i, j, k),
                            grid->Ey(i, j, k),
                            grid->Ez(i, j, k),
                            grid->Bx(i, j, k),
                            grid->By(i, j, k),
                            grid->Bz(i, j, k));

                        func("x"_a = coords.x, "y"_a = coords.y, "z"_a = coords.z,
                            "field_value"_a = std::reference_wrapper<ValueField>(field));

                        grid->Ex(i, j, k) = field.E.x;
                        grid->Ey(i, j, k) = field.E.y;
                        grid->Ez(i, j, k) = field.E.z;

                        grid->Bx(i, j, k) = field.B.x;
                        grid->By(i, j, k) = field.B.y;
                        grid->Bz(i, j, k) = field.B.z;
                    }
        }

        void applyFunction(FunctionPointer _func)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            void(*fValueField)(FP, FP, FP, FP*) = (void(*)(FP, FP, FP, FP*))_func;
            const int chunkSize = 32;
            const int nChunks = grid->numCells.z / chunkSize;
            const int chunkRem = grid->numCells.z % chunkSize;
            const int nx = grid->numCells.x, ny = grid->numCells.y;

            OMP_FOR_COLLAPSE()
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 coords[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
                        FP3 startPosition = grid->ExPosition(i, j, chunk * chunkSize);
                        OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            FP3 position(startPosition.x, startPosition.y,
                                startPosition.z + k * grid->steps.z);
                            coords[k] = derived->convertCoords(position);
                        }
                        OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            int zIndex = k + chunk * chunkSize;
                            ValueField field(grid->Ex(i, j, zIndex),
                                grid->Ey(i, j, zIndex),
                                grid->Ez(i, j, zIndex),
                                grid->Bx(i, j, zIndex),
                                grid->By(i, j, zIndex),
                                grid->Bz(i, j, zIndex));

                            fValueField(coords[k].x, coords[k].y, coords[k].z, &(field.E.x));

                            grid->Ex(i, j, zIndex) = field.E.x;
                            grid->Ey(i, j, zIndex) = field.E.y;
                            grid->Ez(i, j, zIndex) = field.E.z;

                            grid->Bx(i, j, zIndex) = field.B.x;
                            grid->By(i, j, zIndex) = field.B.y;
                            grid->Bz(i, j, zIndex) = field.B.z;
                        }
                    }
        }
    };

    
    // Interface for field solvers with computational grid (grid setters and getters)
    template<class TFieldSolver, class TPyField>
    class pyGridFieldSolverInterface :
        public pySpatialStraggeredFieldInterface<TFieldSolver, TPyField, isFieldSolverSpatialStaggered<TFieldSolver>::value>
    {
    public:

        using TGrid = typename TFieldSolver::GridType;

        void pySetExyz(py::function fEx, py::function fEy, py::function fEz)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            for (int i = 0; i < grid->numCells.x; i++)
                for (int j = 0; j < grid->numCells.y; j++)
                    for (int k = 0; k < grid->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(grid->ExPosition(i, j, k));
                        cEy = derived->convertCoords(grid->EyPosition(i, j, k));
                        cEz = derived->convertCoords(grid->EzPosition(i, j, k));
                        grid->Ex(i, j, k) = fEx("x"_a = cEx.x, "y"_a = cEx.y, "z"_a = cEx.z).template cast<FP>();
                        grid->Ey(i, j, k) = fEy("x"_a = cEy.x, "y"_a = cEy.y, "z"_a = cEy.z).template cast<FP>();
                        grid->Ez(i, j, k) = fEz("x"_a = cEz.x, "y"_a = cEz.y, "z"_a = cEz.z).template cast<FP>();
                    }
        }

        void pySetE(py::function fE)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            for (int i = 0; i < grid->numCells.x; i++)
                for (int j = 0; j < grid->numCells.y; j++)
                    for (int k = 0; k < grid->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(grid->ExPosition(i, j, k));
                        cEy = derived->convertCoords(grid->EyPosition(i, j, k));
                        cEz = derived->convertCoords(grid->EzPosition(i, j, k));
                        grid->Ex(i, j, k) = fE("x"_a = cEx.x, "y"_a = cEx.y, "z"_a = cEx.z).template cast<FP3>().x;
                        grid->Ey(i, j, k) = fE("x"_a = cEy.x, "y"_a = cEy.y, "z"_a = cEy.z).template cast<FP3>().y;
                        grid->Ez(i, j, k) = fE("x"_a = cEz.x, "y"_a = cEz.y, "z"_a = cEz.z).template cast<FP3>().z;
                    }
        }

        void setExyz(FunctionPointer _fEx, FunctionPointer _fEy, FunctionPointer _fEz)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            FP(*fEx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEx;
            FP(*fEy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEy;
            FP(*fEz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEz;
            OMP_FOR_COLLAPSE()
                for (int i = 0; i < grid->numCells.x; i++)
                    for (int j = 0; j < grid->numCells.y; j++)
                        for (int k = 0; k < grid->numCells.z; k++)
                        {
                            FP3 cEx, cEy, cEz;
                            cEx = derived->convertCoords(grid->ExPosition(i, j, k));
                            cEy = derived->convertCoords(grid->EyPosition(i, j, k));
                            cEz = derived->convertCoords(grid->EzPosition(i, j, k));
                            grid->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z);
                            grid->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z);
                            grid->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z);
                        }
        }

        void setExyzt(FunctionPointer _fEx, FunctionPointer _fEy, FunctionPointer _fEz, FP t)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            FP(*fEx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEx;
            FP(*fEy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEy;
            FP(*fEz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEz;
            OMP_FOR_COLLAPSE()
                for (int i = 0; i < grid->numCells.x; i++)
                    for (int j = 0; j < grid->numCells.y; j++)
                        for (int k = 0; k < grid->numCells.z; k++)
                        {
                            FP3 cEx, cEy, cEz;
                            cEx = derived->convertCoords(grid->ExPosition(i, j, k));
                            cEy = derived->convertCoords(grid->EyPosition(i, j, k));
                            cEz = derived->convertCoords(grid->EzPosition(i, j, k));
                            grid->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z, t);
                            grid->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z, t);
                            grid->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z, t);
                        }
        }

        void setE(FunctionPointer _fE)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            FP3(*fE)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fE;
            OMP_FOR_COLLAPSE()
                for (int i = 0; i < grid->numCells.x; i++)
                    for (int j = 0; j < grid->numCells.y; j++)
                        for (int k = 0; k < grid->numCells.z; k++)
                        {
                            FP3 cEx, cEy, cEz;
                            cEx = derived->convertCoords(grid->ExPosition(i, j, k));
                            cEy = derived->convertCoords(grid->EyPosition(i, j, k));
                            cEz = derived->convertCoords(grid->EzPosition(i, j, k));
                            grid->Ex(i, j, k) = fE(cEx.x, cEx.y, cEx.z).x;
                            grid->Ey(i, j, k) = fE(cEy.x, cEy.y, cEy.z).y;
                            grid->Ez(i, j, k) = fE(cEz.x, cEz.y, cEz.z).z;
                        }
        }

        void pySetBxyz(py::function fBx, py::function fBy, py::function fBz)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            for (int i = 0; i < grid->numCells.x; i++)
                for (int j = 0; j < grid->numCells.y; j++)
                    for (int k = 0; k < grid->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(grid->BxPosition(i, j, k));
                        cBy = derived->convertCoords(grid->ByPosition(i, j, k));
                        cBz = derived->convertCoords(grid->BzPosition(i, j, k));
                        grid->Bx(i, j, k) = fBx("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP>();
                        grid->By(i, j, k) = fBy("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP>();
                        grid->Bz(i, j, k) = fBz("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP>();
                    }
        }

        void pySetB(py::function fB)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            for (int i = 0; i < grid->numCells.x; i++)
                for (int j = 0; j < grid->numCells.y; j++)
                    for (int k = 0; k < grid->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(grid->BxPosition(i, j, k));
                        cBy = derived->convertCoords(grid->ByPosition(i, j, k));
                        cBz = derived->convertCoords(grid->BzPosition(i, j, k));
                        grid->Bx(i, j, k) = fB("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP3>().x;
                        grid->By(i, j, k) = fB("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP3>().y;
                        grid->Bz(i, j, k) = fB("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP3>().z;
                    }
        }

        void setBxyz(FunctionPointer _fBx, FunctionPointer _fBy, FunctionPointer _fBz)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            FP(*fBx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBx;
            FP(*fBy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBy;
            FP(*fBz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBz;
            OMP_FOR_COLLAPSE()
                for (int i = 0; i < grid->numCells.x; i++)
                    for (int j = 0; j < grid->numCells.y; j++)
                        for (int k = 0; k < grid->numCells.z; k++)
                        {
                            FP3 cBx, cBy, cBz;
                            cBx = derived->convertCoords(grid->BxPosition(i, j, k));
                            cBy = derived->convertCoords(grid->ByPosition(i, j, k));
                            cBz = derived->convertCoords(grid->BzPosition(i, j, k));
                            grid->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z);
                            grid->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z);
                            grid->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z);
                        }
        }

        void setBxyzt(FunctionPointer _fBx, FunctionPointer _fBy, FunctionPointer _fBz, FP t)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            FP(*fBx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBx;
            FP(*fBy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBy;
            FP(*fBz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBz;
            OMP_FOR_COLLAPSE()
                for (int i = 0; i < grid->numCells.x; i++)
                    for (int j = 0; j < grid->numCells.y; j++)
                        for (int k = 0; k < grid->numCells.z; k++)
                        {
                            FP3 cBx, cBy, cBz;
                            cBx = derived->convertCoords(grid->BxPosition(i, j, k));
                            cBy = derived->convertCoords(grid->ByPosition(i, j, k));
                            cBz = derived->convertCoords(grid->BzPosition(i, j, k));
                            grid->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z, t);
                            grid->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z, t);
                            grid->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z, t);
                        }
        }

        void setB(FunctionPointer _fB)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            FP3(*fB)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fB;
            OMP_FOR_COLLAPSE()
                for (int i = 0; i < grid->numCells.x; i++)
                    for (int j = 0; j < grid->numCells.y; j++)
                        for (int k = 0; k < grid->numCells.z; k++)
                        {
                            FP3 cBx, cBy, cBz;
                            cBx = derived->convertCoords(grid->BxPosition(i, j, k));
                            cBy = derived->convertCoords(grid->ByPosition(i, j, k));
                            cBz = derived->convertCoords(grid->BzPosition(i, j, k));
                            grid->Bx(i, j, k) = fB(cBx.x, cBx.y, cBx.z).x;
                            grid->By(i, j, k) = fB(cBy.x, cBy.y, cBy.z).y;
                            grid->Bz(i, j, k) = fB(cBz.x, cBz.y, cBz.z).z;
                        }
        }

        void pySetJxyz(py::function fJx, py::function fJy, py::function fJz)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            for (int i = 0; i < grid->numCells.x; i++)
                for (int j = 0; j < grid->numCells.y; j++)
                    for (int k = 0; k < grid->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(grid->JxPosition(i, j, k));
                        cJy = derived->convertCoords(grid->JyPosition(i, j, k));
                        cJz = derived->convertCoords(grid->JzPosition(i, j, k));
                        grid->Jx(i, j, k) = fJx("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP>();
                        grid->Jy(i, j, k) = fJy("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP>();
                        grid->Jz(i, j, k) = fJz("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP>();
                    }
        }

        void pySetJ(py::function fJ)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            for (int i = 0; i < grid->numCells.x; i++)
                for (int j = 0; j < grid->numCells.y; j++)
                    for (int k = 0; k < grid->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(grid->JxPosition(i, j, k));
                        cJy = derived->convertCoords(grid->JyPosition(i, j, k));
                        cJz = derived->convertCoords(grid->JzPosition(i, j, k));
                        grid->Jx(i, j, k) = fJ("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP3>().x;
                        grid->Jy(i, j, k) = fJ("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP3>().y;
                        grid->Jz(i, j, k) = fJ("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP3>().z;
                    }
        }

        void setJxyz(FunctionPointer _fJx, FunctionPointer _fJy, FunctionPointer _fJz)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            FP(*fJx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJx;
            FP(*fJy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJy;
            FP(*fJz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJz;
            OMP_FOR_COLLAPSE()
                for (int i = 0; i < grid->numCells.x; i++)
                    for (int j = 0; j < grid->numCells.y; j++)
                        for (int k = 0; k < grid->numCells.z; k++)
                        {
                            FP3 cJx, cJy, cJz;
                            cJx = derived->convertCoords(grid->JxPosition(i, j, k));
                            cJy = derived->convertCoords(grid->JyPosition(i, j, k));
                            cJz = derived->convertCoords(grid->JzPosition(i, j, k));
                            grid->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z);
                            grid->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z);
                            grid->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z);
                        }
        }

        void setJxyzt(FunctionPointer _fJx, FunctionPointer _fJy, FunctionPointer _fJz, FP t)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            FP(*fJx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJx;
            FP(*fJy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJy;
            FP(*fJz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJz;
            OMP_FOR_COLLAPSE()
                for (int i = 0; i < grid->numCells.x; i++)
                    for (int j = 0; j < grid->numCells.y; j++)
                        for (int k = 0; k < grid->numCells.z; k++)
                        {
                            FP3 cJx, cJy, cJz;
                            cJx = derived->convertCoords(grid->JxPosition(i, j, k));
                            cJy = derived->convertCoords(grid->JyPosition(i, j, k));
                            cJz = derived->convertCoords(grid->JzPosition(i, j, k));
                            grid->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z, t);
                            grid->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z, t);
                            grid->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z, t);
                        }
        }

        void setJ(FunctionPointer _fJ)
        {
            TPyField* derived = static_cast<TPyField*>(this);
            TGrid* grid = derived->getField()->getGrid();
            TFieldSolver* fieldSolver = derived->getField()->getFieldSolver();
            FP3(*fJ)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fJ;
            OMP_FOR_COLLAPSE()
                for (int i = 0; i < grid->numCells.x; i++)
                    for (int j = 0; j < grid->numCells.y; j++)
                        for (int k = 0; k < grid->numCells.z; k++)
                        {
                            FP3 cJx, cJy, cJz;
                            cJx = derived->convertCoords(grid->JxPosition(i, j, k));
                            cJy = derived->convertCoords(grid->JyPosition(i, j, k));
                            cJz = derived->convertCoords(grid->JzPosition(i, j, k));
                            grid->Jx(i, j, k) = fJ(cJx.x, cJx.y, cJx.z).x;
                            grid->Jy(i, j, k) = fJ(cJy.x, cJy.y, cJy.z).y;
                            grid->Jz(i, j, k) = fJ(cJz.x, cJz.y, cJz.z).z;
                        }
        }

        FP3 getE(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getE(coords);
        }
        FP3 getB(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getB(coords);
        }
        FP3 getJ(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getJ(coords);
        }

        FP getEx(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getEx(coords);
        }
        FP getEy(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getEy(coords);
        }
        FP getEz(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getEz(coords);
        }

        FP getBx(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getBx(coords);
        }
        FP getBy(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getBy(coords);
        }
        FP getBz(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getBz(coords);
        }

        FP getJx(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getJx(coords);
        }
        FP getJy(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getJy(coords);
        }
        FP getJz(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getJz(coords);
        }

        void getFields(const FP3& coords, FP3& e, FP3& b) const {
            static_cast<const TPyField*>(this)->getField()->getGrid()->getFields(coords, e, b);
        }
    };


    // Interface for solvers with Poisson's equation or not
    template <class TFieldSolver, class TPyField, bool>
    class pyPoissonFieldSolverInterface {};

    // Interface for solvers that support Poisson's equation
    template <class TFieldSolver, class TPyField>
    class pyPoissonFieldSolverInterface<TFieldSolver, TPyField, true> {
    public:
        void convertFieldsPoissonEquation() {
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->convertFieldsPoissonEquation();
        }
    };


    // Interface for solvers with field generator or not
    template <class TFieldSolver, class TPyField, bool>
    class pyFieldGeneratorSolverInterface {};
    
    // Interface for solvers that support field generator
    template <class TFieldSolver, class TPyField>
    class pyFieldGeneratorSolverInterface<TFieldSolver, TPyField, true> {
    public:

        using FunctionType = typename TFieldSolver::FieldGeneratorType::FunctionType;
        
        void setFieldGenerator(
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FunctionPointer bxFunc, FunctionPointer byFunc, FunctionPointer bzFunc,
            FunctionPointer exFunc, FunctionPointer eyFunc, FunctionPointer ezFunc,
            bool isXLeftBorderEnabled, bool isYLeftBorderEnabled, bool isZLeftBorderEnabled,
            bool isXRightBorderEnabled, bool isYRightBorderEnabled, bool isZRightBorderEnabled
        )
        {
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->setFieldGenerator(
                leftGenIndex, rightGenIndex,
                cppFunc(bxFunc), cppFunc(byFunc), cppFunc(bzFunc),
                cppFunc(exFunc), cppFunc(eyFunc), cppFunc(ezFunc),
                Int3(isXLeftBorderEnabled, isYLeftBorderEnabled, isZLeftBorderEnabled),
                Int3(isXRightBorderEnabled, isYRightBorderEnabled, isZRightBorderEnabled)
            );
        }

        void setFieldGeneratorAllFunctions(
            const Int3& leftGenIndex, const Int3& rightGenIndex,

            FunctionPointer xMinBx, FunctionPointer xMinBy, FunctionPointer xMinBz,
            FunctionPointer xMaxBx, FunctionPointer xMaxBy, FunctionPointer xMaxBz,
            FunctionPointer yMinBx, FunctionPointer yMinBy, FunctionPointer yMinBz,
            FunctionPointer yMaxBx, FunctionPointer yMaxBy, FunctionPointer yMaxBz,
            FunctionPointer zMinBx, FunctionPointer zMinBy, FunctionPointer zMinBz,
            FunctionPointer zMaxBx, FunctionPointer zMaxBy, FunctionPointer zMaxBz,

            FunctionPointer xMinEx, FunctionPointer xMinEy, FunctionPointer xMinEz,
            FunctionPointer xMaxEx, FunctionPointer xMaxEy, FunctionPointer xMaxEz,
            FunctionPointer yMinEx, FunctionPointer yMinEy, FunctionPointer yMinEz,
            FunctionPointer yMaxEx, FunctionPointer yMaxEy, FunctionPointer yMaxEz,
            FunctionPointer zMinEx, FunctionPointer zMinEy, FunctionPointer zMinEz,
            FunctionPointer zMaxEx, FunctionPointer zMaxEy, FunctionPointer zMaxEz,

            bool isXLeftBorderEnabled, bool isYLeftBorderEnabled, bool isZLeftBorderEnabled,
            bool isXRightBorderEnabled, bool isYRightBorderEnabled, bool isZRightBorderEnabled
        )
        {
            std::array<std::array<FunctionType, 3>, 3> leftBFunc;
            leftBFunc[0][0] = cppFunc(xMinBx); leftBFunc[0][1] = cppFunc(xMinBy); leftBFunc[0][2] = cppFunc(xMinBz);
            leftBFunc[1][0] = cppFunc(yMinBx); leftBFunc[1][1] = cppFunc(yMinBy); leftBFunc[1][2] = cppFunc(yMinBz);
            leftBFunc[2][0] = cppFunc(zMinBx); leftBFunc[2][1] = cppFunc(zMinBy); leftBFunc[2][2] = cppFunc(zMinBz);

            std::array<std::array<FunctionType, 3>, 3> rightBFunc;
            rightBFunc[0][0] = cppFunc(xMaxBx); rightBFunc[0][1] = cppFunc(xMaxBy); rightBFunc[0][2] = cppFunc(xMaxBz);
            rightBFunc[1][0] = cppFunc(yMaxBx); rightBFunc[1][1] = cppFunc(yMaxBy); rightBFunc[1][2] = cppFunc(yMaxBz);
            rightBFunc[2][0] = cppFunc(zMaxBx); rightBFunc[2][1] = cppFunc(zMaxBy); rightBFunc[2][2] = cppFunc(zMaxBz);

            std::array<std::array<FunctionType, 3>, 3> leftEFunc;
            leftEFunc[0][0] = cppFunc(xMinEx); leftEFunc[0][1] = cppFunc(xMinEy); leftEFunc[0][2] = cppFunc(xMinEz);
            leftEFunc[1][0] = cppFunc(yMinEx); leftEFunc[1][1] = cppFunc(yMinEy); leftEFunc[1][2] = cppFunc(yMinEz);
            leftEFunc[2][0] = cppFunc(zMinEx); leftEFunc[2][1] = cppFunc(zMinEy); leftEFunc[2][2] = cppFunc(zMinEz);

            std::array<std::array<FunctionType, 3>, 3> rightEFunc;
            rightEFunc[0][0] = cppFunc(xMaxEx); rightEFunc[0][1] = cppFunc(xMaxEy); rightEFunc[0][2] = cppFunc(xMaxEz);
            rightEFunc[1][0] = cppFunc(yMaxEx); rightEFunc[1][1] = cppFunc(yMaxEy); rightEFunc[1][2] = cppFunc(yMaxEz);
            rightEFunc[2][0] = cppFunc(zMaxEx); rightEFunc[2][1] = cppFunc(zMaxEy); rightEFunc[2][2] = cppFunc(zMaxEz);

            static_cast<TPyField*>(this)->getField()->getFieldSolver()->setFieldGenerator(
                leftGenIndex, rightGenIndex,
                leftBFunc, rightBFunc, leftEFunc, rightEFunc,
                Int3(isXLeftBorderEnabled, isYLeftBorderEnabled, isZLeftBorderEnabled),
                Int3(isXRightBorderEnabled, isYRightBorderEnabled, isZRightBorderEnabled)
            );
        }

    protected:

        FunctionType cppFunc(FunctionPointer f) {
            return FunctionType((FP(*)(FP, FP, FP, FP))f);
        };

    };


    // Interface for solvers with periodical boundary conditions or not
    template <class TFieldSolver, class TPyField, bool>
    class pyPeriodicalBoundaryConditionSolverInterface {};

    // Interface for solvers that support periodical boundary conditions
    template <class TFieldSolver, class TPyField>
    class pyPeriodicalBoundaryConditionSolverInterface<TFieldSolver, TPyField, true> {
    public:

        using BoundaryConditionType = typename TFieldSolver::PeriodicalBoundaryConditionType;

        void setPeriodicalBoundaryCondition() {
            static_cast<TPyField*>(this)->getField()
                ->getFieldSolver()->template setBoundaryCondition<BoundaryConditionType>();
        }

        void setPeriodicalBoundaryCondition(CoordinateEnum axis) {
            static_cast<TPyField*>(this)->getField()
                ->getFieldSolver()->template setBoundaryCondition<BoundaryConditionType>(axis);
        }
    };


    // Interface for solvers with reflect boundary conditions or not
    template <class TFieldSolver, class TPyField, bool>
    class pyReflectBoundaryConditionSolverInterface {};

    // Interface for solvers that support reflect boundary conditions
    template <class TFieldSolver, class TPyField>
    class pyReflectBoundaryConditionSolverInterface<TFieldSolver, TPyField, true> {
    public:

        using BoundaryConditionType = typename TFieldSolver::ReflectBoundaryConditionType;

        void setReflectBoundaryCondition() {
            static_cast<TPyField*>(this)->getField()
                ->getFieldSolver()->template setBoundaryCondition<BoundaryConditionType>();
        }

        void setReflectBoundaryCondition(CoordinateEnum axis) {
            static_cast<TPyField*>(this)->getField()
                ->getFieldSolver()->template setBoundaryCondition<BoundaryConditionType>(axis);
        }
    };


    // Interface for solvers with pml or not
    template <class TFieldSolver, class TPyField, bool>
    class pyPMLSolverInterface {};

    // Interface for solvers that support PML
    template <class TFieldSolver, class TPyField>
    class pyPMLSolverInterface<TFieldSolver, TPyField, true> {
    public:
        void setPML(int sizePMLx, int sizePMLy, int sizePMLz) {
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->setPML(sizePMLx, sizePMLy, sizePMLz);
        }
    };


    // Interface for analytical and non-analytical fields
    template<class TFieldSolver, class TPyField, bool ifAnalyticalField>
    class pyAnalyticalFieldInterface {};

    // Interface for numerical fields
    template<class TFieldSolver, class TPyField>
    class pyAnalyticalFieldInterface<TFieldSolver, TPyField, false> :
        public pyGridFieldSolverInterface<TFieldSolver, TPyField>,
        public pyPoissonFieldSolverInterface<TFieldSolver, TPyField, isFieldSolverSupportPoissonEquation<TFieldSolver>::value>,
        public pyFieldGeneratorSolverInterface<TFieldSolver, TPyField, isFieldSolverSupportFieldGenerator<TFieldSolver>::value>,
        public pyPeriodicalBoundaryConditionSolverInterface<TFieldSolver, TPyField, isFieldSolverSupportPeriodicalBoundaryCondition<TFieldSolver>::value>,
        public pyReflectBoundaryConditionSolverInterface<TFieldSolver, TPyField, isFieldSolverSupportReflectBoundaryCondition<TFieldSolver>::value>,
        public pyPMLSolverInterface<TFieldSolver, TPyField, isFieldSolverSupportPml<TFieldSolver>::value>
    {};
    
    // Interface for analytical fields
    template<class TFieldSolver, class TPyField>
    class pyAnalyticalFieldInterface<TFieldSolver, TPyField, true> {
    public:

        void setExyz(FunctionPointer fEx, FunctionPointer fEy, FunctionPointer fEz) {
            FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEx;
            FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEy;
            FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEz;
            static_cast<const TPyField*>(this)->getField()->getGrid()->setE(fx, fy, fz);
        }

        void setBxyz(FunctionPointer fBx, FunctionPointer fBy, FunctionPointer fBz) {
            FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBx;
            FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBy;
            FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBz;
            static_cast<const TPyField*>(this)->getField()->getGrid()->setB(fx, fy, fz);
        }

        void setJxyz(FunctionPointer fJx, FunctionPointer fJy, FunctionPointer fJz) {
            FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fJx;
            FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fJy;
            FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fJz;
            static_cast<const TPyField*>(this)->getField()->getGrid()->setJ(fx, fy, fz);
        }

        FP3 getE(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getE(coords);
        }
        FP3 getB(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getB(coords);
        }
        FP3 getJ(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getJ(coords);
        }

        FP getEx(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getEx(coords);
        }
        FP getEy(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getEy(coords);
        }
        FP getEz(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getEz(coords);
        }

        FP getBx(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getBx(coords);
        }
        FP getBy(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getBy(coords);
        }
        FP getBz(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getBz(coords);
        }

        FP getJx(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getJx(coords);
        }
        FP getJy(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getJy(coords);
        }
        FP getJz(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getJz(coords);
        }

        FP3 getEt(FP x, FP y, FP z, FP t) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getE(x, y, z, t);
        }
        FP3 getBt(FP x, FP y, FP z, FP t) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getB(x, y, z, t);
        }
        FP3 getJt(FP x, FP y, FP z, FP t) const {
            return static_cast<const TPyField*>(this)->getField()->getGrid()->getJ(x, y, z, t);
        }
    };


    // Common interface for all fields
    template<class TFieldSolver, class TPyField>
    class pyFieldInterface :
        public pyAnalyticalFieldInterface<TFieldSolver, TPyField, isFieldAnalytical<TFieldSolver>::value>
    {
    public:

        void setTime(FP time) {
            static_cast<TPyField*>(this)->getField()->setTime(time);
        }

        FP getTime() {
            return static_cast<TPyField*>(this)->getField()->getTime();
        }

        void changeTimeStep(double dt) {
            static_cast<TPyField*>(this)->getField()->setTimeStep(dt);
        }

        void updateFields() {
            static_cast<TPyField*>(this)->getField()->updateFields();
        }

        void advance(FP dt) {
            static_cast<TPyField*>(this)->getField()->advance(dt);
        }

        void refresh() {
            static_cast<TPyField*>(this)->getField()->refresh();
        }
    };

}
