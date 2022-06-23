#pragma once
#include <memory>

#include "Grid.h"
#include "AnalyticalField.h"
#include "FieldValue.h"
#include "Mapping.h"
#include "Fdtd.h"
#include "Psatd.h"
#include "Pstd.h"
#include "Mapping.h"
#include "Field.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace pfc
{

    template <class TGrid, class TFieldSolver, class TPyField, bool ifStraggered>
    class pyStraggeredFieldInterface {};

    // Interface for spatial straggered grids
    template <class TGrid, class TFieldSolver, class TPyField>
    class pyStraggeredFieldInterface<TGrid, TFieldSolver, TPyField, true>
    {
    public:

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
                            cEx[k] = derived->convertCoords(grid->ExPosition(i, j, chunk * chunkSize),
                                fieldSolver->timeShiftE);
                            cEy[k] = derived->convertCoords(grid->EyPosition(i, j, chunk * chunkSize),
                                fieldSolver->timeShiftE);
                            cEz[k] = derived->convertCoords(grid->EzPosition(i, j, chunk * chunkSize),
                                fieldSolver->timeShiftE);

                            cBx[k] = derived->convertCoords(grid->BxPosition(i, j, chunk * chunkSize),
                                fieldSolver->timeShiftB);
                            cBy[k] = derived->convertCoords(grid->ByPosition(i, j, chunk * chunkSize),
                                fieldSolver->timeShiftB);
                            cBz[k] = derived->convertCoords(grid->BzPosition(i, j, chunk * chunkSize),
                                fieldSolver->timeShiftB);
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
    template <class TGrid, class TFieldSolver, class TPyField>
    class pyStraggeredFieldInterface<TGrid, TFieldSolver, TPyField, false>
    {
    public:

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

        void setEMField(int64_t _fValueField)
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

        void applyFunction(int64_t _func)
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


    template<class TGrid, class TFieldSolver, class TPyField, bool ifAnalyticalField>
    class pyGridFieldInterface {};

    // Interface for fields with computational grid
    template<class TGrid, class TFieldSolver, class TPyField>
    class pyGridFieldInterface<TGrid, TFieldSolver, TPyField, false> :
        public pyStraggeredFieldInterface<TGrid, TFieldSolver,
        TPyField, TGrid::ifFieldsSpatialStraggered && TGrid::ifFieldsTimeStraggered>
    {
    public:

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
                        cEx = derived->convertCoords(grid->ExPosition(i, j, k), fieldSolver->timeShiftE);
                        cEy = derived->convertCoords(grid->EyPosition(i, j, k), fieldSolver->timeShiftE);
                        cEz = derived->convertCoords(grid->EzPosition(i, j, k), fieldSolver->timeShiftE);
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
                        cEx = derived->convertCoords(grid->ExPosition(i, j, k), fieldSolver->timeShiftE);
                        cEy = derived->convertCoords(grid->EyPosition(i, j, k), fieldSolver->timeShiftE);
                        cEz = derived->convertCoords(grid->EzPosition(i, j, k), fieldSolver->timeShiftE);
                        grid->Ex(i, j, k) = fE("x"_a = cEx.x, "y"_a = cEx.y, "z"_a = cEx.z).template cast<FP3>().x;
                        grid->Ey(i, j, k) = fE("x"_a = cEy.x, "y"_a = cEy.y, "z"_a = cEy.z).template cast<FP3>().y;
                        grid->Ez(i, j, k) = fE("x"_a = cEz.x, "y"_a = cEz.y, "z"_a = cEz.z).template cast<FP3>().z;
                    }
        }

        void setExyz(int64_t _fEx, int64_t _fEy, int64_t _fEz)
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
                        cEx = derived->convertCoords(grid->ExPosition(i, j, k), fieldSolver->timeShiftE);
                        cEy = derived->convertCoords(grid->EyPosition(i, j, k), fieldSolver->timeShiftE);
                        cEz = derived->convertCoords(grid->EzPosition(i, j, k), fieldSolver->timeShiftE);
                        grid->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z);
                        grid->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z);
                        grid->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z);
                    }
        }

        void setExyzt(int64_t _fEx, int64_t _fEy, int64_t _fEz, FP t)
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
                        cEx = derived->convertCoords(grid->ExPosition(i, j, k), fieldSolver->timeShiftE);
                        cEy = derived->convertCoords(grid->EyPosition(i, j, k), fieldSolver->timeShiftE);
                        cEz = derived->convertCoords(grid->EzPosition(i, j, k), fieldSolver->timeShiftE);
                        grid->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z, t + fieldSolver->timeShiftE);
                        grid->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z, t + fieldSolver->timeShiftE);
                        grid->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z, t + fieldSolver->timeShiftE);
                    }
        }

        void setE(int64_t _fE)
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
                        cEx = derived->convertCoords(grid->ExPosition(i, j, k), fieldSolver->timeShiftE);
                        cEy = derived->convertCoords(grid->EyPosition(i, j, k), fieldSolver->timeShiftE);
                        cEz = derived->convertCoords(grid->EzPosition(i, j, k), fieldSolver->timeShiftE);
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
                        cBx = derived->convertCoords(grid->BxPosition(i, j, k), fieldSolver->timeShiftB);
                        cBy = derived->convertCoords(grid->ByPosition(i, j, k), fieldSolver->timeShiftB);
                        cBz = derived->convertCoords(grid->BzPosition(i, j, k), fieldSolver->timeShiftB);
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
                        cBx = derived->convertCoords(grid->BxPosition(i, j, k), fieldSolver->timeShiftB);
                        cBy = derived->convertCoords(grid->ByPosition(i, j, k), fieldSolver->timeShiftB);
                        cBz = derived->convertCoords(grid->BzPosition(i, j, k), fieldSolver->timeShiftB);
                        grid->Bx(i, j, k) = fB("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP3>().x;
                        grid->By(i, j, k) = fB("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP3>().y;
                        grid->Bz(i, j, k) = fB("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP3>().z;
                    }
        }

        void setBxyz(int64_t _fBx, int64_t _fBy, int64_t _fBz)
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
                        cBx = derived->convertCoords(grid->BxPosition(i, j, k), fieldSolver->timeShiftB);
                        cBy = derived->convertCoords(grid->ByPosition(i, j, k), fieldSolver->timeShiftB);
                        cBz = derived->convertCoords(grid->BzPosition(i, j, k), fieldSolver->timeShiftB);
                        grid->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z);
                        grid->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z);
                        grid->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z);
                    }
        }

        void setBxyzt(int64_t _fBx, int64_t _fBy, int64_t _fBz, FP t)
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
                        cBx = derived->convertCoords(grid->BxPosition(i, j, k), fieldSolver->timeShiftB);
                        cBy = derived->convertCoords(grid->ByPosition(i, j, k), fieldSolver->timeShiftB);
                        cBz = derived->convertCoords(grid->BzPosition(i, j, k), fieldSolver->timeShiftB);
                        grid->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z, t + fieldSolver->timeShiftB);
                        grid->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z, t + fieldSolver->timeShiftB);
                        grid->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z, t + fieldSolver->timeShiftB);
                    }
        }

        void setB(int64_t _fB)
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
                        cBx = derived->convertCoords(grid->BxPosition(i, j, k), fieldSolver->timeShiftB);
                        cBy = derived->convertCoords(grid->ByPosition(i, j, k), fieldSolver->timeShiftB);
                        cBz = derived->convertCoords(grid->BzPosition(i, j, k), fieldSolver->timeShiftB);
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
                        cJx = derived->convertCoords(grid->JxPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJy = derived->convertCoords(grid->JyPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJz = derived->convertCoords(grid->JzPosition(i, j, k), fieldSolver->timeShiftJ);
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
                        cJx = derived->convertCoords(grid->JxPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJy = derived->convertCoords(grid->JyPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJz = derived->convertCoords(grid->JzPosition(i, j, k), fieldSolver->timeShiftJ);
                        grid->Jx(i, j, k) = fJ("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP3>().x;
                        grid->Jy(i, j, k) = fJ("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP3>().y;
                        grid->Jz(i, j, k) = fJ("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP3>().z;
                    }
        }

        void setJxyz(int64_t _fJx, int64_t _fJy, int64_t _fJz)
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
                        cJx = derived->convertCoords(grid->JxPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJy = derived->convertCoords(grid->JyPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJz = derived->convertCoords(grid->JzPosition(i, j, k), fieldSolver->timeShiftJ);
                        grid->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z);
                        grid->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z);
                        grid->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z);
                    }
        }

        void setJxyzt(int64_t _fJx, int64_t _fJy, int64_t _fJz, FP t)
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
                        cJx = derived->convertCoords(grid->JxPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJy = derived->convertCoords(grid->JyPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJz = derived->convertCoords(grid->JzPosition(i, j, k), fieldSolver->timeShiftJ);
                        grid->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z, t + fieldSolver->timeShiftJ);
                        grid->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z, t + fieldSolver->timeShiftJ);
                        grid->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z, t + fieldSolver->timeShiftJ);
                    }
        }

        void setJ(int64_t _fJ)
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
                        cJx = derived->convertCoords(grid->JxPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJy = derived->convertCoords(grid->JyPosition(i, j, k), fieldSolver->timeShiftJ);
                        cJz = derived->convertCoords(grid->JzPosition(i, j, k), fieldSolver->timeShiftJ);
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


    // Interface for analytical fields
    template<class TGrid, class TFieldSolver, class TPyField>
    class pyGridFieldInterface<TGrid, TFieldSolver, TPyField, true> {
    public:

        void setExyz(int64_t fEx, int64_t fEy, int64_t fEz) {
            FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEx;
            FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEy;
            FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEz;
            static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->setE(fx, fy, fz);
        }

        void setBxyz(int64_t fBx, int64_t fBy, int64_t fBz) {
            FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBx;
            FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBy;
            FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBz;
            static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->setB(fx, fy, fz);
        }

        void setJxyz(int64_t fJx, int64_t fJy, int64_t fJz) {
            FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fJx;
            FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fJy;
            FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fJz;
            static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->setJ(fx, fy, fz);
        }

        FP3 getE(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getE(coords);
        }
        FP3 getB(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getB(coords);
        }
        FP3 getJ(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getJ(coords);
        }

        FP getEx(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getEx(coords);
        }
        FP getEy(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getEy(coords);
        }
        FP getEz(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getEz(coords);
        }

        FP getBx(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getBx(coords);
        }
        FP getBy(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getBy(coords);
        }
        FP getBz(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getBz(coords);
        }

        FP getJx(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getJx(coords);
        }
        FP getJy(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getJy(coords);
        }
        FP getJz(const FP3& coords) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getJz(coords);
        }

        FP3 getEt(FP x, FP y, FP z, FP t) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getE(x, y, z, t);
        }
        FP3 getBt(FP x, FP y, FP z, FP t) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getB(x, y, z, t);
        }
        FP3 getJt(FP x, FP y, FP z, FP t) const {
            return static_cast<const TPyField*>(this)->getField()->getAnalyticalField()->getJ(x, y, z, t);
        }
    };


    template <class TGrid, class TFieldSolver, class TPyField, bool>
    class pyPoissonFieldSolverInterface {};

    // Interface for solvers that support Poisson's equation
    template <class TGrid, class TFieldSolver, class TPyField>
    class pyPoissonFieldSolverInterface<TGrid, TFieldSolver, TPyField, true> {
    public:
        void convertFieldsPoissonEquation() {
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->convertFieldsPoissonEquation();
        }
    };


    template <class TGrid, class TFieldSolver, class TPyField, bool>
    class pyFieldGeneratorSolverInterface {};

    // Interface for solvers that support field generator
    template <class TGrid, class TFieldSolver, class TPyField>
    class pyFieldGeneratorSolverInterface<TGrid, TFieldSolver, TPyField, true> {
    public:
        void setPeriodicalFieldGenerator() {
            PeriodicalFieldGenerator<TGrid::gridType> gen;
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->setFieldGenerator(&gen);
        }

        void setReflectFieldGenerator() {
            ReflectFieldGenerator<TGrid::gridType> gen;
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->setFieldGenerator(&gen);
        }
    };


    template <class TGrid, class TFieldSolver, class TPyField, bool>
    class pyPMLSolverInterface {};

    // Interface for solvers that support PML
    template <class TGrid, class TFieldSolver, class TPyField>
    class pyPMLSolverInterface<TGrid, TFieldSolver, TPyField, true> {
    public:
        void setPML(int sizePMLx, int sizePMLy, int sizePMLz) {
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->setPML(sizePMLx, sizePMLy, sizePMLz);
        }
    };


    // Common interface for field solvers
    template<class TGrid, class TFieldSolver, class TPyField>
    class pyFieldSolverInterface :
        public pyPoissonFieldSolverInterface<TGrid, TFieldSolver, TPyField,
        std::is_same<TFieldSolver, PSATD>::value || std::is_same<TFieldSolver, PSATDPoisson>::value ||
        std::is_same<TFieldSolver, PSATDTimeStraggered>::value ||
        std::is_same<TFieldSolver, PSATDTimeStraggeredPoisson>::value>,
        public pyFieldGeneratorSolverInterface<TGrid, TFieldSolver, TPyField,
        std::is_same<TFieldSolver, FDTD>::value>,
        public pyPMLSolverInterface<TGrid, TFieldSolver, TPyField,
        !std::is_same<TFieldSolver, NoFieldSolver>::value>
    {
    public:

        void setTime(FP time) {
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->globalTime = time;
        }

        FP getTime() {
            return static_cast<TPyField*>(this)->getField()->getFieldSolver()->globalTime;
        }

        void changeTimeStep(double dt) {
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->setTimeStep(dt);
        }

        void updateFields() {
            static_cast<TPyField*>(this)->getField()->getFieldSolver()->updateFields();
        }

        void advance(FP dt) {
            auto fieldSolver = static_cast<TPyField*>(this)->getField()->getFieldSolver();
            FP oldDt = fieldSolver->dt;
            fieldSolver->setTimeStep(dt);
            fieldSolver->updateFields();
            fieldSolver->setTimeStep(oldDt);
        }
    };


    // Common interface for all field
    template<class TGrid, class TFieldSolver, class TPyField>
    class pyFieldInterface :
        public pyGridFieldInterface<TGrid, TFieldSolver, TPyField,
        std::is_same<TFieldSolver, NoFieldSolver>::value>,
        public pyFieldSolverInterface<TGrid, TFieldSolver, TPyField>
    {
    public:

        void refresh() {
            static_cast<TPyField*>(this)->getField()->refresh();
        }
    };

}
