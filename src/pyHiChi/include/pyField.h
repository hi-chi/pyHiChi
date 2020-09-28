#pragma once
#include <memory>
#include "Grid.h"
#include "FieldValue.h"
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
    class pyStraggeredFieldIntarface {};

    // spatial straggered grids
    template <class TGrid, class TFieldSolver, class TDerived>
    class pyStraggeredFieldIntarface<TGrid, TFieldSolver, TDerived, true>
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
    class pyStraggeredFieldIntarface<TGrid, TFieldSolver, TDerived, false>
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
    class pyFieldGridInterface : public pyStraggeredFieldIntarface<TGrid, TFieldSolver,
        TDerived, TGrid::ifFieldsSpatialStraggered && TGrid::ifFieldsTimeStraggered>
    {
    public:

        pyFieldGridInterface()
        {
            fEt[0] = 0; fEt[1] = 0; fEt[2] = 0;
            fBt[0] = 0; fBt[1] = 0; fBt[2] = 0;
            isAnalytical = false;
        }

        void setAnalytical(int64_t _fEx, int64_t _fEy, int64_t _fEz, int64_t _fBx, int64_t _fBy, int64_t _fBz)
        {
            fEt[0] = _fEx; fEt[1] = _fEy; fEt[2] = _fEz;
            fBt[0] = _fBx; fBt[1] = _fBy; fBt[2] = _fBz;
            isAnalytical = true;
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

        FP3 getE(const FP3& coords) const
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

        FP3 getB(const FP3& coords) const
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

        FP3 getJ(const FP3& coords) const {
            return static_cast<const TDerived*>(this)->getFieldEntity()->getJ(coords);
        }

        void getFields(const FP3& coords, FP3& e, FP3& b) const {
            static_cast<const TDerived*>(this)->getFieldEntity()->getFields(coords, e, b);
        }

    private:

        int64_t fEt[3], fBt[3];
        bool isAnalytical;

    };


    template <class TGrid, class TFieldSolver, class TDerived, bool>
    class pyPoissonFieldSolverInterface {};

    template <class TGrid, class TFieldSolver, class TDerived>
    class pyPoissonFieldSolverInterface<TGrid, TFieldSolver, TDerived, true> {
    public:
        void convertFieldsPoissonEquation() {
            static_cast<TDerived*>(this)->getFieldEntity()->convertFieldsPoissonEquation();
        }
    };

    template <class TGrid, class TFieldSolver, class TDerived>
    class pyPoissonFieldSolverInterface<TGrid, TFieldSolver, TDerived, false> {
    public:
        void convertFieldsPoissonEquation() {
            std::cout
                << "WARNING: the used field does not include the 'convertFieldsPoissonEquation' method"
                << std::endl;
        }
    };


    template <class TGrid, class TFieldSolver, class TDerived, bool>
    class pyFieldGeneratorSolverInterface {};

    template <class TGrid, class TFieldSolver, class TDerived>
    class pyFieldGeneratorSolverInterface<TGrid, TFieldSolver, TDerived, true> {
    public:
        void setFieldGenerator(FieldGenerator<TGrid::gridType>* generator) {
            static_cast<TDerived*>(this)->getFieldEntity()->setFieldGenerator(generator);
        }
    };

    template <class TGrid, class TFieldSolver, class TDerived>
    class pyFieldGeneratorSolverInterface<TGrid, TFieldSolver, TDerived, false> {
    public:
        void setFieldGenerator(FieldGenerator<TGrid::gridType>* generator) {
            std::cout
                << "WARNING: the used field does not include the 'setFieldGenerator' method"
                << std::endl;
        }
    };


    template<class TGrid, class TFieldSolver, class TDerived>
    class pyFieldSolverInterface :
        public pyPoissonFieldSolverInterface<TGrid, TFieldSolver, TDerived,
        std::is_same<TFieldSolver, PSATD>::value || std::is_same<TFieldSolver, PSATDPoisson>::value ||
        std::is_same<TFieldSolver, PSATDTimeStraggered>::value ||
        std::is_same<TFieldSolver, PSATDTimeStraggeredPoisson>::value>,
        public pyFieldGeneratorSolverInterface<TGrid, TFieldSolver, TDerived,
        std::is_same<TFieldSolver, FDTD>::value>
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


    template<class TGrid, class TFieldSolver, class TDerived>
    class pyFieldInterface:
        public pyFieldGridInterface<TGrid, TFieldSolver, TDerived>,
        public pyFieldSolverInterface<TGrid, TFieldSolver, TDerived>
    {
    public:

        using BaseGridInterface =
            pyFieldGridInterface<TGrid, TFieldSolver, TDerived>;
        using BaseSolverInterface =
            pyFieldSolverInterface<TGrid, TFieldSolver, TDerived>;

        TGrid* getGrid() const {
            return static_cast<TGrid*>(static_cast<const TDerived*>(this)->getFieldEntity());
        }

        TFieldSolver* getFieldSolver() const {
            return static_cast<TFieldSolver*>(static_cast<const TDerived*>(this)->getFieldEntity());
        }

        void refresh() {
            static_cast<TDerived*>(this)->getFieldEntity()->refresh();
        }
    };


    class pyFieldBase {
    public:

        virtual FP3 getE(const FP3& coords) const = 0;
        virtual FP3 getB(const FP3& coords) const = 0;
        virtual FP3 getJ(const FP3& coords) const = 0;

        void getFields(const FP3& coords, FP3& e, FP3& b) const {
            e = getE(coords);
            b = getB(coords);
        }

        virtual void updateFields() = 0;
        virtual void advance(FP dt) = 0;

        virtual std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const = 0;
    };


    template<class TGrid, class TFieldSolver>
    class pyField : public pyFieldInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>,
        public pyFieldBase
    {
        using BaseInterface = pyFieldInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>;

    public:

        pyField(const Int3 & numInternalCells,
            const FP3 & minCoords, const FP3 & steps, FP dt) :
            fieldEntity(new pyFieldEntity<TGrid, TFieldSolver>(numInternalCells,
                minCoords, steps, dt))
        {}

        pyField(const std::shared_ptr<pyField<TGrid, TFieldSolver>>& other,
            const std::shared_ptr<Mapping>& mapping) :
            pyWrappedField(other), mapping(mapping)
        {}

        inline pyFieldEntity<TGrid, TFieldSolver>* getFieldEntity() const {
            if (fieldEntity) return fieldEntity.get();
            return pyWrappedField->getFieldEntity();
        }

        inline FP3 convertCoords(const FP3& coords, FP timeShift = 0.0) const {
            bool status = true;
            return getDirectCoords(coords, getFieldEntity()->globalTime + timeShift, &status);
        }

        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pyField<TGrid, TFieldSolver>>(
                    std::static_pointer_cast<pyField<TGrid, TFieldSolver>>(self), mapping
                    )
                );
        }

        inline FP3 getE(const FP3& coords) const override {
            bool status = true;
            FP time = getFieldEntity()->globalTime + getFieldEntity()->timeShiftE;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0, 0, 0);
            return BaseInterface::getE(inverseCoords);
        }

        inline FP3 getB(const FP3& coords) const override {
            bool status = true;
            FP time = getFieldEntity()->globalTime + getFieldEntity()->timeShiftB;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0, 0, 0);
            return BaseInterface::getB(inverseCoords);
        }

        FP3 getJ(const FP3& coords) const override {
            bool status = true;
            FP time = getFieldEntity()->globalTime + getFieldEntity()->timeShiftJ;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0, 0, 0);
            return BaseInterface::getJ(inverseCoords);
        }

        void updateFields() override {
            return BaseInterface::updateFields();
        }

        void advance(FP dt) override {
            return BaseInterface::advance(dt);
        }

    protected:

        inline FP3 getDirectCoords(const FP3& coords, FP time, bool* status) const {
            FP3 coords_ = coords;
            *status = true;
            if (pyWrappedField) coords_ = pyWrappedField->getDirectCoords(coords_, time, status);
            bool status2 = true;
            if (mapping) coords_ = mapping->getDirectCoords(coords_, time, &status2);
            *status = *status && status2;
            return coords_;
        }

        inline FP3 getInverseCoords(const FP3& coords, FP time, bool* status) const {
            FP3 coords_ = coords;
            *status = true;
            if (pyWrappedField) coords_ = pyWrappedField->getInverseCoords(coords_, time, status);
            bool status2 = true;
            if (mapping) coords_ = mapping->getInverseCoords(coords_, time, &status2);
            *status = *status && status2;
            return coords_;
        }

    private:

        // the simple grid state
        // if fieldEntity!=0 then pyField is a memory owner
        std::unique_ptr<pyFieldEntity<TGrid, TFieldSolver>> fieldEntity;

        // the mapping grid state
        // if pyWrappedField!=0 then pyField is a wrapper
        std::shared_ptr<pyField<TGrid, TFieldSolver>> pyWrappedField;
        std::shared_ptr<Mapping> mapping;

    };

    typedef pyField<YeeGrid, FDTD> pyYeeField;
    typedef pyField<PSTDGrid, PSTD> pyPSTDField;
    typedef pyField<PSATDGrid, PSATD> pyPSATDField;
    typedef pyField<PSATDGrid, PSATDPoisson> pyPSATDPoissonField;
    typedef pyField<PSATDTimeStraggeredGrid, PSATDTimeStraggered> pyPSATDTimeStraggeredField;
    typedef pyField<PSATDTimeStraggeredGrid, PSATDTimeStraggeredPoisson> pyPSATDTimeStraggeredPoissonField;


    class pySumField : public pyFieldBase
    {
    public:

        pySumField(const std::shared_ptr<pyFieldBase>& pyWrappedField1,
            const std::shared_ptr<pyFieldBase>& pyWrappedField2) :
            pyWrappedField1(pyWrappedField1), pyWrappedField2(pyWrappedField2)
        {}

        pySumField(const std::shared_ptr<pySumField>& other,
            const std::shared_ptr<Mapping>& mapping) :
            pyWrappedField1(other->pyWrappedField1->applyMapping(other->pyWrappedField1, mapping)),
            pyWrappedField2(other->pyWrappedField2->applyMapping(other->pyWrappedField2, mapping))
        {}

        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pySumField>(
                    std::static_pointer_cast<pySumField>(self), mapping
                    )
                );
        }

        FP3 getE(const FP3& coords) const override {
            return pyWrappedField1->getE(coords) + pyWrappedField2->getE(coords);
        }

        FP3 getB(const FP3& coords) const override {
            return pyWrappedField1->getB(coords) + pyWrappedField2->getB(coords);
        }

        FP3 getJ(const FP3& coords) const override {
            return pyWrappedField1->getJ(coords) + pyWrappedField2->getJ(coords);
        }

        void updateFields() override {
            pyWrappedField1->updateFields();
            pyWrappedField2->updateFields();
        }

        void advance(FP dt) override {
            pyWrappedField1->advance(dt);
            pyWrappedField2->advance(dt);
        }

    private:

        std::shared_ptr<pyFieldBase> pyWrappedField1;
        std::shared_ptr<pyFieldBase> pyWrappedField2;
    };


    class pyMulField : public pyFieldBase {
    public:

        pyMulField(const std::shared_ptr<pyFieldBase>& pyWrappedField, FP factor) :
            pyWrappedField(pyWrappedField),
            factor(factor)
        {}

        pyMulField(const std::shared_ptr<pyMulField>& other,
            const std::shared_ptr<Mapping>& mapping) :
            pyWrappedField(other->pyWrappedField->applyMapping(other->pyWrappedField, mapping))
        {}

        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pyMulField>(
                    std::static_pointer_cast<pyMulField>(self), mapping
                    )
                );
        }

        FP3 getE(const FP3& coords) const override {
            return pyWrappedField->getE(coords) * factor;
        }

        FP3 getB(const FP3& coords) const override {
            return pyWrappedField->getB(coords) * factor;
        }

        FP3 getJ(const FP3& coords) const override {
            return pyWrappedField->getJ(coords) * factor;
        }

        void updateFields() override {
            pyWrappedField->updateFields();
        }

        void advance(FP dt) override {
            pyWrappedField->advance(dt);
        }

    private:

        FP factor = 1.0;
        std::shared_ptr<pyFieldBase> pyWrappedField;
    };

}