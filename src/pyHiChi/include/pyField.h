#pragma once
#include <memory>
#include "Grid.h"
#include "Mapping.h"
#include "Fdtd.h"
#include "Psatd.h"
#include "Pstd.h"
#include "Mapping.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace pfc
{
    template <class TGrid, class TFieldSolver>
    class pyFieldEntity : public TGrid, public TFieldSolver {
    public:

        pyFieldEntity(const Int3 & numInternalCells,
            const FP3 & minCoords, const FP3 & steps, FP dt) :
            TGrid(Int3(numInternalCells), minCoords, steps, numInternalCells),
            TFieldSolver(static_cast<TGrid*>(this), dt)
        {}

        void refresh() {
            this->globalTime = 0.0;
        }
    };

    template <class TGrid, class TFieldSolver, class TDerived, bool ifStraggered>
    class pyStraggeredFieldSetterIntarface {};

    // spatial straggered grids
    template <class TGrid, class TFieldSolver, class TDerived>
    class pyStraggeredFieldSetterIntarface<TGrid, TFieldSolver, TDerived, true>
    {
    public:

        template <class FieldConfigurationType>
        void setFieldConfiguration(const FieldConfigurationType* fieldConf) {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            const int chunkSize = 32;
            const int nChunks = fieldEntity->numCells.z / chunkSize;
            const int chunkRem = fieldEntity->numCells.z % chunkSize;
            const int nx = fieldEntity->numCells.x, ny = fieldEntity->numCells.y;

#pragma omp parallel for collapse(2)
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 cEx[chunkSize], cEy[chunkSize], cEz[chunkSize];
                        FP3 cBx[chunkSize], cBy[chunkSize], cBz[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
#pragma ivdep
                        for (int k = 0; k < kLast; k++) {
                            cEx[k] = derived->convertCoords(fieldEntity->ExPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftE);
                            cEy[k] = derived->convertCoords(fieldEntity->EyPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftE);
                            cEz[k] = derived->convertCoords(fieldEntity->EzPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftE);

                            cBx[k] = derived->convertCoords(fieldEntity->BxPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftB);
                            cBy[k] = derived->convertCoords(fieldEntity->ByPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftB);
                            cBz[k] = derived->convertCoords(fieldEntity->BzPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftB);
                        }
#pragma ivdep
#pragma omp simd
                        for (int k = 0; k < kLast; k++) {
                            fieldEntity->Ex(i, j, k) = fieldConf->getE(cEx[k].x, cEx[k].y, cEx[k].z).x;
                            fieldEntity->Ey(i, j, k) = fieldConf->getE(cEy[k].x, cEy[k].y, cEy[k].z).y;
                            fieldEntity->Ez(i, j, k) = fieldConf->getE(cEz[k].x, cEz[k].y, cEz[k].z).z;

                            fieldEntity->Bx(i, j, k) = fieldConf->getB(cBx[k].x, cBx[k].y, cBx[k].z).x;
                            fieldEntity->By(i, j, k) = fieldConf->getB(cBy[k].x, cBy[k].y, cBy[k].z).y;
                            fieldEntity->Bz(i, j, k) = fieldConf->getB(cBz[k].x, cBz[k].y, cBz[k].z).z;
                        }
                    }
        }

    };


    // collocated grids
    template <class TGrid, class TFieldSolver, class TDerived>
    class pyStraggeredFieldSetterIntarface<TGrid, TFieldSolver, TDerived, false>
    {
    public:

        template <class FieldConfigurationType>
        void setFieldConfiguration(const FieldConfigurationType* fieldConf) {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            const int chunkSize = 32;
            const int nChunks = fieldEntity->numCells.z / chunkSize;
            const int chunkRem = fieldEntity->numCells.z % chunkSize;
            const int nx = fieldEntity->numCells.x, ny = fieldEntity->numCells.y;

#pragma omp parallel for collapse(2)
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 coords[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
                        FP3 startPosition = fieldEntity->ExPosition(i, j, chunk * chunkSize);
#pragma ivdep
                        for (int k = 0; k < kLast; k++) {
                            FP3 position(startPosition.x, startPosition.y, startPosition.z + k * fieldEntity->steps.z);
                            coords[k] = derived->convertCoords(position);
                        }
#pragma ivdep
#pragma omp simd
                        for (int k = 0; k < kLast; k++) {
                            FP3 E, B;
                            fieldConf->getEB(coords[k].x, coords[k].y, coords[k].z, &E, &B);

                            fieldEntity->Ex(i, j, k + chunk * chunkSize) = E.x;
                            fieldEntity->Ey(i, j, k + chunk * chunkSize) = E.y;
                            fieldEntity->Ez(i, j, k + chunk * chunkSize) = E.z;

                            fieldEntity->Bx(i, j, k + chunk * chunkSize) = B.x;
                            fieldEntity->By(i, j, k + chunk * chunkSize) = B.y;
                            fieldEntity->Bz(i, j, k + chunk * chunkSize) = B.z;
                        }
                    }
        }

        void pySetEMField(py::function fValueField)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 coords = derived->convertCoords(fieldEntity->ExPosition(i, j, k));
                        ValueField field = fValueField("x"_a = coords.x, "y"_a = coords.y, "z"_a = coords.z).
                            template cast<ValueField>();

                        fieldEntity->Ex(i, j, k) = field.E.x;
                        fieldEntity->Ey(i, j, k) = field.E.y;
                        fieldEntity->Ez(i, j, k) = field.E.z;

                        fieldEntity->Bx(i, j, k) = field.B.x;
                        fieldEntity->By(i, j, k) = field.B.y;
                        fieldEntity->Bz(i, j, k) = field.B.z;
                    }
        }

        void setEMField(int64_t _fValueField)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            void(*fValueField)(FP, FP, FP, FP*) = (void(*)(FP, FP, FP, FP*))_fValueField;
            const int chunkSize = 32;
            const int nChunks = fieldEntity->numCells.z / chunkSize;
            const int chunkRem = fieldEntity->numCells.z % chunkSize;
            const int nx = fieldEntity->numCells.x, ny = fieldEntity->numCells.y;

#pragma omp parallel for collapse(2)
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 coords[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
                        FP3 startPosition = fieldEntity->ExPosition(i, j, chunk * chunkSize);
#pragma ivdep
                        for (int k = 0; k < kLast; k++) {
                            FP3 position(startPosition.x, startPosition.y,
                                startPosition.z + k * fieldEntity->steps.z);
                            coords[k] = derived->convertCoords(position);
                        }
#pragma ivdep
#pragma omp simd
                        for (int k = 0; k < kLast; k++) {
                            ValueField field(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                            fValueField(coords[k].x, coords[k].y, coords[k].z, &(field.E.x));

                            fieldEntity->Ex(i, j, k + chunk * chunkSize) = field.E.x;
                            fieldEntity->Ey(i, j, k + chunk * chunkSize) = field.E.y;
                            fieldEntity->Ez(i, j, k + chunk * chunkSize) = field.E.z;

                            fieldEntity->Bx(i, j, k + chunk * chunkSize) = field.B.x;
                            fieldEntity->By(i, j, k + chunk * chunkSize) = field.B.y;
                            fieldEntity->Bz(i, j, k + chunk * chunkSize) = field.B.z;
                        }
                    }
        }
    };


    template<class TGrid, class TFieldSolver, class TDerived>
    class pyFieldSetterInterface : public pyStraggeredFieldSetterIntarface<TGrid, TFieldSolver,
        TDerived, TGrid::ifFieldsSpatialStraggered && TGrid::ifFieldsTimeStraggered>
    {
    public:

        pyFieldSetterInterface()
        {
            fEt[0] = 0; fEt[1] = 0; fEt[2] = 0;
            fBt[0] = 0; fBt[1] = 0; fBt[2] = 0;
            isAnalytical = false;
        }

        void analyticalUpdateFields(FP t)
        {
            if (isAnalytical)
            {
                setExyzt(fEt[0], fEt[1], fEt[2], t);
                setBxyzt(fBt[0], fBt[1], fBt[2], t);
            }
        }

        void pySetExyz(py::function fEx, py::function fEy, py::function fEz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(fieldEntity->ExPosition(i, j, k), fieldEntity->timeShiftE);
                        cEy = derived->convertCoords(fieldEntity->EyPosition(i, j, k), fieldEntity->timeShiftE);
                        cEz = derived->convertCoords(fieldEntity->EzPosition(i, j, k), fieldEntity->timeShiftE);
                        fieldEntity->Ex(i, j, k) = fEx("x"_a = cEx.x, "y"_a = cEx.y, "z"_a = cEx.z).template cast<FP>();
                        fieldEntity->Ey(i, j, k) = fEy("x"_a = cEy.x, "y"_a = cEy.y, "z"_a = cEy.z).template cast<FP>();
                        fieldEntity->Ez(i, j, k) = fEz("x"_a = cEz.x, "y"_a = cEz.y, "z"_a = cEz.z).template cast<FP>();
                    }
        }

        void pySetE(py::function fE)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(fieldEntity->ExPosition(i, j, k), fieldEntity->timeShiftE);
                        cEy = derived->convertCoords(fieldEntity->EyPosition(i, j, k), fieldEntity->timeShiftE);
                        cEz = derived->convertCoords(fieldEntity->EzPosition(i, j, k), fieldEntity->timeShiftE);
                        fieldEntity->Ex(i, j, k) = fE("x"_a = cEx.x, "y"_a = cEx.y, "z"_a = cEx.z).template cast<FP3>().x;
                        fieldEntity->Ey(i, j, k) = fE("x"_a = cEy.x, "y"_a = cEy.y, "z"_a = cEy.z).template cast<FP3>().y;
                        fieldEntity->Ez(i, j, k) = fE("x"_a = cEz.x, "y"_a = cEz.y, "z"_a = cEz.z).template cast<FP3>().z;
                    }
        }

        void setExyz(int64_t _fEx, int64_t _fEy, int64_t _fEz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fEx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEx;
            FP(*fEy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEy;
            FP(*fEz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEz;
#pragma omp parallel for
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(fieldEntity->ExPosition(i, j, k), fieldEntity->timeShiftE);
                        cEy = derived->convertCoords(fieldEntity->EyPosition(i, j, k), fieldEntity->timeShiftE);
                        cEz = derived->convertCoords(fieldEntity->EzPosition(i, j, k), fieldEntity->timeShiftE);
                        fieldEntity->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z);
                        fieldEntity->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z);
                        fieldEntity->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z);
                    }
        }

        void setExyzt(int64_t _fEx, int64_t _fEy, int64_t _fEz, FP t)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fEx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEx;
            FP(*fEy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEy;
            FP(*fEz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEz;
#pragma omp parallel for
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(fieldEntity->ExPosition(i, j, k), fieldEntity->timeShiftE);
                        cEy = derived->convertCoords(fieldEntity->EyPosition(i, j, k), fieldEntity->timeShiftE);
                        cEz = derived->convertCoords(fieldEntity->EzPosition(i, j, k), fieldEntity->timeShiftE);
                        fieldEntity->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z, t + fieldEntity->timeShiftE);
                        fieldEntity->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z, t + fieldEntity->timeShiftE);
                        fieldEntity->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z, t + fieldEntity->timeShiftE);
                    }
        }

        void setE(int64_t _fE)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP3(*fE)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fE;
#pragma omp parallel for
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(fieldEntity->ExPosition(i, j, k), fieldEntity->timeShiftE);
                        cEy = derived->convertCoords(fieldEntity->EyPosition(i, j, k), fieldEntity->timeShiftE);
                        cEz = derived->convertCoords(fieldEntity->EzPosition(i, j, k), fieldEntity->timeShiftE);
                        fieldEntity->Ex(i, j, k) = fE(cEx.x, cEx.y, cEx.z).x;
                        fieldEntity->Ey(i, j, k) = fE(cEy.x, cEy.y, cEy.z).y;
                        fieldEntity->Ez(i, j, k) = fE(cEz.x, cEz.y, cEz.z).z;
                    }
        }

        void pySetBxyz(py::function fBx, py::function fBy, py::function fBz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(fieldEntity->BxPosition(i, j, k), fieldEntity->timeShiftB);
                        cBy = derived->convertCoords(fieldEntity->ByPosition(i, j, k), fieldEntity->timeShiftB);
                        cBz = derived->convertCoords(fieldEntity->BzPosition(i, j, k), fieldEntity->timeShiftB);
                        fieldEntity->Bx(i, j, k) = fBx("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP>();
                        fieldEntity->By(i, j, k) = fBy("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP>();
                        fieldEntity->Bz(i, j, k) = fBz("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP>();
                    }
        }

        void pySetB(py::function fB)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(fieldEntity->BxPosition(i, j, k), fieldEntity->timeShiftB);
                        cBy = derived->convertCoords(fieldEntity->ByPosition(i, j, k), fieldEntity->timeShiftB);
                        cBz = derived->convertCoords(fieldEntity->BzPosition(i, j, k), fieldEntity->timeShiftB);
                        fieldEntity->Bx(i, j, k) = fB("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP3>().x;
                        fieldEntity->By(i, j, k) = fB("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP3>().y;
                        fieldEntity->Bz(i, j, k) = fB("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP3>().z;
                    }
        }

        void setBxyz(int64_t _fBx, int64_t _fBy, int64_t _fBz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fBx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBx;
            FP(*fBy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBy;
            FP(*fBz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBz;
#pragma omp parallel for
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(fieldEntity->BxPosition(i, j, k), fieldEntity->timeShiftB);
                        cBy = derived->convertCoords(fieldEntity->ByPosition(i, j, k), fieldEntity->timeShiftB);
                        cBz = derived->convertCoords(fieldEntity->BzPosition(i, j, k), fieldEntity->timeShiftB);
                        fieldEntity->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z);
                        fieldEntity->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z);
                        fieldEntity->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z);
                    }
        }

        void setBxyzt(int64_t _fBx, int64_t _fBy, int64_t _fBz, FP t)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fBx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBx;
            FP(*fBy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBy;
            FP(*fBz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBz;
#pragma omp parallel for
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(fieldEntity->BxPosition(i, j, k), fieldEntity->timeShiftB);
                        cBy = derived->convertCoords(fieldEntity->ByPosition(i, j, k), fieldEntity->timeShiftB);
                        cBz = derived->convertCoords(fieldEntity->BzPosition(i, j, k), fieldEntity->timeShiftB);
                        fieldEntity->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z, t + fieldEntity->timeShiftB);
                        fieldEntity->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z, t + fieldEntity->timeShiftB);
                        fieldEntity->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z, t + fieldEntity->timeShiftB);
                    }
        }

        void setB(int64_t _fB)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP3(*fB)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fB;
#pragma omp parallel for
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(fieldEntity->BxPosition(i, j, k), fieldEntity->timeShiftB);
                        cBy = derived->convertCoords(fieldEntity->ByPosition(i, j, k), fieldEntity->timeShiftB);
                        cBz = derived->convertCoords(fieldEntity->BzPosition(i, j, k), fieldEntity->timeShiftB);
                        fieldEntity->Bx(i, j, k) = fB(cBx.x, cBx.y, cBx.z).x;
                        fieldEntity->By(i, j, k) = fB(cBy.x, cBy.y, cBy.z).y;
                        fieldEntity->Bz(i, j, k) = fB(cBz.x, cBz.y, cBz.z).z;
                    }
        }

        void pySetJxyz(py::function fJx, py::function fJy, py::function fJz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJx("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP>();
                        fieldEntity->Jy(i, j, k) = fJy("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP>();
                        fieldEntity->Jz(i, j, k) = fJz("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP>();
                    }
        }

        void pySetJ(py::function fJ)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJ("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP3>().x;
                        fieldEntity->Jy(i, j, k) = fJ("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP3>().y;
                        fieldEntity->Jz(i, j, k) = fJ("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP3>().z;
                    }
        }

        void setJxyz(int64_t _fJx, int64_t _fJy, int64_t _fJz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fJx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJx;
            FP(*fJy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJy;
            FP(*fJz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJz;
#pragma omp parallel for
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z);
                        fieldEntity->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z);
                        fieldEntity->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z);
                    }
        }

        void setJxyzt(int64_t _fJx, int64_t _fJy, int64_t _fJz, FP t)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fJx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJx;
            FP(*fJy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJy;
            FP(*fJz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJz;
#pragma omp parallel for
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z, t + fieldEntity->timeShiftJ);
                        fieldEntity->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z, t + fieldEntity->timeShiftJ);
                        fieldEntity->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z, t + fieldEntity->timeShiftJ);
                    }
        }

        void setJ(int64_t _fJ)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP3(*fJ)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fJ;
#pragma omp parallel for
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJ(cJx.x, cJx.y, cJx.z).x;
                        fieldEntity->Jy(i, j, k) = fJ(cJy.x, cJy.y, cJy.z).y;
                        fieldEntity->Jz(i, j, k) = fJ(cJz.x, cJz.y, cJz.z).z;
                    }
        }

    protected:

        void setAnalytical(int64_t _fEx, int64_t _fEy, int64_t _fEz, int64_t _fBx, int64_t _fBy, int64_t _fBz)
        {
            fEt[0] = _fEx; fEt[1] = _fEy; fEt[2] = _fEz;
            fBt[0] = _fBx; fBt[1] = _fBy; fBt[2] = _fBz;
            isAnalytical = true;
        }

    private:

        int64_t fEt[3], fBt[3];
        bool isAnalytical;

    };


    template<class TGrid, class TFieldSolver, class TDerived>
    class pyFieldGetterInterface
    {
    public:

        pyFieldGetterInterface()
        {
            fEt[0] = 0; fEt[1] = 0; fEt[2] = 0;
            fBt[0] = 0; fBt[1] = 0; fBt[2] = 0;
            isAnalytical = false;
        }

        FP3 getE(const FP3& coords) const {
            return static_cast<TDerived*>(this)->getE(coords);
        }

        FP3 getB(const FP3& coords) const {
            return static_cast<TDerived*>(this)->getB(coords);
        }

        FP3 getJ(const FP3& coords) const {
            return static_cast<TDerived*>(this)->getJ(coords);
        }

        void getFields(const FP3& coords, FP3 & e, FP3 & b) const {
            static_cast<const TDerived*>(this)->getFieldEntity()->getFields(coords, e, b);
        }

    protected:

        FP3 getEByDirectCoords(const FP3& coords) const
        {

            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity =
                static_cast<const TDerived*>(this)->getFieldEntity();
            FP3 result;
            if (isAnalytical)
            {
                FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEt[0];
                FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEt[1];
                FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEt[2];
                FP time = fieldEntity->globalTime + fieldEntity->timeShiftE;
                result[0] = fx(coords.x, coords.y, coords.z, time);
                result[1] = fy(coords.x, coords.y, coords.z, time);
                result[2] = fz(coords.x, coords.y, coords.z, time);
            }
            else {
                result = fieldEntity->getE(coords);
            }
            return result;
        }

        FP3 getBByDirectCoords(const FP3& coords) const
        {
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity =
                static_cast<const TDerived*>(this)->getFieldEntity();
            FP3 result;
            if (isAnalytical)
            {
                FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBt[0];
                FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBt[1];
                FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBt[2];
                FP time = fieldEntity->globalTime + fieldEntity->timeShiftB;
                result[0] = fx(coords.x, coords.y, coords.z, time);
                result[1] = fy(coords.x, coords.y, coords.z, time);
                result[2] = fz(coords.x, coords.y, coords.z, time);
            }
            else {
                result = fieldEntity->getB(coords);
            }
            return result;
        }

        FP3 getJByDirectCoords(const FP3& coords) const {
            return static_cast<const TDerived*>(this)->getFieldEntity()->getJ(coords);
        }

        void setAnalytical(int64_t _fEx, int64_t _fEy, int64_t _fEz, int64_t _fBx, int64_t _fBy, int64_t _fBz)
        {
            fEt[0] = _fEx; fEt[1] = _fEy; fEt[2] = _fEz;
            fBt[0] = _fBx; fBt[1] = _fBy; fBt[2] = _fBz;
            isAnalytical = true;
        }

    private:

        int64_t fEt[3], fBt[3];
        bool isAnalytical;

    };


    template<class TGrid, class TFieldSolver, class TDerived>
    class pyFieldSolverInterface
    {
    public:

        void setTime(FP time) {
            static_cast<TDerived*>(this)->getFieldEntity()->globalTime = time;
        }

        FP getTime() {
            return static_cast<TDerived*>(this)->getFieldEntity()->globalTime;
        }

        void setPML(int sizePMLx, int sizePMLy, int sizePMLz) {
            static_cast<TDerived*>(this)->getFieldEntity()->setPML(sizePMLx, sizePMLy, sizePMLz);
        }

        void setFieldGenerator(FieldGenerator<TGrid::gridType>* generator) {
            static_cast<TDerived*>(this)->getFieldEntity()->setFieldGenerator(generator);
        }

        void changeTimeStep(double dt) {
            static_cast<TDerived*>(this)->getFieldEntity()->setTimeStep(dt);
        }

        void updateFields() {
            static_cast<TDerived*>(this)->getFieldEntity()->updateFields();
        }

        void advance(FP dt) {
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = static_cast<TDerived*>(this)->getFieldEntity();
            FP oldDt = fieldEntity->dt;
            fieldEntity->setTimeStep(dt);
            fieldEntity->updateFields();
            fieldEntity->setTimeStep(oldDt);
        }
    };


    template<class TGrid, class TFieldSolver>
    class pyField :
        public pyFieldSetterInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>,
        public pyFieldGetterInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>,
        public pyFieldSolverInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>
    {
        using BaseSetterInterface = pyFieldSetterInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>;
        using BaseGetterInterface = pyFieldGetterInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>;
        using BaseSolverInterface = pyFieldSolverInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>;

    public:

        pyField(const Int3 & numInternalCells,
            const FP3 & minCoords, const FP3 & steps, FP dt) :
            fieldEntity(std::make_unique<pyFieldEntity<TGrid, TFieldSolver>>(numInternalCells,
                minCoords, steps, dt))
        {}

        pyField(const std::shared_ptr<pyField<TGrid, TFieldSolver>>& pyFieldRef,
            const std::shared_ptr<Mapping>& mapping) :
            pyFieldRef(pyFieldRef), mapping(mapping)
        {}

        void setAnalytical(int64_t fEx, int64_t fEy, int64_t fEz, int64_t fBx, int64_t fBy, int64_t fBz)
        {
            BaseSetterInterface::setAnalytical(fEx, fEy, fEz, fBx, fBy, fBz);
            BaseGetterInterface::setAnalytical(fEx, fEy, fEz, fBx, fBy, fBz);
        }

        static std::shared_ptr<pyField<TGrid, TFieldSolver>> applyMapping(
            const std::shared_ptr<pyField<TGrid, TFieldSolver>>& pyFieldRef,
            const std::shared_ptr<Mapping>& mapping) {
            return std::make_shared<pyField<TGrid, TFieldSolver>>(pyFieldRef, mapping);
        }

        inline pyFieldEntity<TGrid, TFieldSolver>* getFieldEntity() const {
            if (fieldEntity) return fieldEntity.get();
            return pyFieldRef->getFieldEntity();
        }

        TGrid* getGrid() const {
            return static_cast<TGrid*>(getFieldEntity());
        }

        TFieldSolver* getFieldSolver() const {
            return static_cast<TFieldSolver*>(getFieldEntity());
        }

        inline FP3 convertCoords(const FP3& coords, FP timeShift = 0.0) const {
            bool status = true;
            return getDirectCoords(coords, getFieldEntity()->globalTime + timeShift, &status);
        }

        inline FP3 getE(const FP3& coords) const {
            bool status = true;
            FP time = getFieldEntity()->globalTime + getFieldEntity()->timeShiftE;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0, 0, 0);
            return BaseGetterInterface::getEByDirectCoords(inverseCoords);
        }

        inline FP3 getB(const FP3& coords) const {
            bool status = true;
            FP time = getFieldEntity()->globalTime + getFieldEntity()->timeShiftB;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0, 0, 0);
            return BaseGetterInterface::getBByDirectCoords(inverseCoords);
        }

        inline FP3 getJ(const FP3& coords) const {
            bool status = true;
            FP time = getFieldEntity()->globalTime + getFieldEntity()->timeShiftJ;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0, 0, 0);
            return BaseGetterInterface::getJByDirectCoords(inverseCoords);
        }

        void refresh() {
            getFieldEntity()->refresh();
        }

    protected:

        inline FP3 getDirectCoords(const FP3& coords, FP time, bool* status) const {
            FP3 coords_ = coords;
            *status = true;
            if (pyFieldRef) coords_ = pyFieldRef->getDirectCoords(coords_, time, status);
            bool status2 = true;
            if (mapping) coords_ = mapping->getDirectCoords(coords_, time, &status2);
            *status = *status && status2;
            return coords_;
        }

        inline FP3 getInverseCoords(const FP3& coords, FP time, bool* status) const {
            FP3 coords_ = coords;
            *status = true;
            if (pyFieldRef) coords_ = pyFieldRef->getInverseCoords(coords_, time, status);
            bool status2 = true;
            if (mapping) coords_ = mapping->getInverseCoords(coords_, time, &status2);
            *status = *status && status2;
            return coords_;
        }

    private:

        // the simple grid state
        std::unique_ptr<pyFieldEntity<TGrid, TFieldSolver>> fieldEntity;  // if !=0 then pyField is a memory owner

        // the mapping grid state
        std::shared_ptr<pyField<TGrid, TFieldSolver>> pyFieldRef;  // if !=0 then pyField is a wrapper
        std::shared_ptr<Mapping> mapping;

    };

    typedef pyField<YeeGrid, FDTD> pyYeeField;
    typedef pyField<PSTDGrid, PSTD> pyPSTDField;
    typedef pyField<PSATDGrid, PSATD> pyPSATDField;
    typedef pyField<PSATDGrid, PSATDPoisson> pyPSATDPoissonField;
    typedef pyField<PSATDTimeStraggeredGrid, PSATDTimeStraggered> pyPSATDTimeStraggeredField;
    typedef pyField<PSATDTimeStraggeredGrid, PSATDTimeStraggeredPoisson> pyPSATDTimeStraggeredPoissonField;

}