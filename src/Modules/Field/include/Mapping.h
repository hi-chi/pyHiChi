#pragma once
#include "Vectors.h"
#include "Constants.h"
#include "Enums.h"
#include <functional>
#include <array>

namespace pfc {

    class Mapping {

    public:

        // coordinate transformations
        // inverse transform
        virtual FP3 getDirectCoords(const FP3& coords, FP time = 0.0, bool* status = 0) {
            setOkStatus(status);
            return coords;
        }
        // direct transform
        virtual FP3 getInverseCoords(const FP3& coords, FP time = 0.0, bool* status = 0) {
            setOkStatus(status);
            return coords;
        }

        // sometimes we need to transform field too, for example, when rotating fields
        virtual bool isRequireFieldTransform() {
            return false;
        }
        // inverse transform
        virtual FP3 getDirectFields(const FP3& fields, FP time = 0.0, bool status = true) {
            if (status) return fields;
            return FP3(0.0, 0.0, 0.0);
        }
        // direct transform
        virtual FP3 getInverseFields(const FP3& fields, FP time = 0.0, bool status = true) {
            if (status) return fields;
            return FP3(0.0, 0.0, 0.0);
        }

        virtual Mapping* createInstance() = 0;

    protected:

        void setFailStatus(bool* status) {
            if (status) *status = false;
        }

        void setOkStatus(bool* status) {
            if (status) *status = true;
        }

    };


    class IdentityMapping : public Mapping {

    public:

        // create identity mapping on a segment [a,b)
        IdentityMapping(FP3 a, FP3 b) : a(a), b(b) {}

        FP3 getDirectCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            if (status) *status = (coords >= a && coords < b) ? true : false;
            return coords;
        }

        FP3 getInverseCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            if (status) *status = (coords >= a && coords < b) ? true : false;
            return coords;
        }

        Mapping* createInstance() override {
            return new IdentityMapping(*this);
        }

        FP3 a, b;

    };


    class PeriodicalMapping : public Mapping {

    public:

        // create periodical mapping: axis = ...[cMin, cMax)[cMin, cMax)[cMin, cMax)...
        PeriodicalMapping(CoordinateEnum axis, FP cMin, FP cMax) :
            axis(axis), cMin(cMin), cMax(cMax), D(cMax-cMin) {}

        FP3 getDirectCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            if (status) *status = (coords[(int)axis] >= cMin && coords[(int)axis] < cMax) ? true : false;
            return coords;
        }

        FP3 getInverseCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            setOkStatus(status);
            FP3 inverseCoords = coords;
            double tmp;
            FP frac = std::modf((coords[(int)axis] - cMin) / D, &tmp);
            inverseCoords[(int)axis] = cMin + (frac >= 0 ? frac : (1 + frac)) * D;
            return inverseCoords;
        }

        Mapping* createInstance() override {
            return new PeriodicalMapping(*this);
        }

        FP cMin, cMax, D;
        CoordinateEnum axis = CoordinateEnum::x;

    };


    class RotationMapping : public Mapping {

    public:

        static const int dim = 3;
        using TRotMatrix = std::array<std::array<FP, (size_t)dim>, (size_t)dim>;

        // rotation around the axis according to the rule of the right hand
        RotationMapping(CoordinateEnum axis, FP angle) {
            createRotationMatrix(rotationMatrix, (int)axis, angle);
        }

        // rotation around the axis + change of polarization vector
        RotationMapping(CoordinateEnum axis, FP angle, CoordinateEnum propDir, FP polarizationAngle) {
            TRotMatrix polarizationMatrix = { 0.0 }, coordRotMatrix = { 0.0 };
            createRotationMatrix(polarizationMatrix, (int)propDir, polarizationAngle);
            createRotationMatrix(coordRotMatrix, (int)axis, angle);
            mulMatrices(coordRotMatrix, polarizationMatrix, this->rotationMatrix);
        }

        FP3 getDirectCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            setOkStatus(status);
            return mulRotationMatrix(coords);
        }

        FP3 getInverseCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            setOkStatus(status);
            return mulInverseRotationMatrix(coords);
        }

        virtual bool isRequireFieldTransform() {
            return true;
        }

        FP3 getDirectFields(const FP3& fields, FP time = 0.0, bool status = true) override {
            if (status) return mulInverseRotationMatrix(fields);
            return FP3(0.0, 0.0, 0.0);
        }

        FP3 getInverseFields(const FP3& fields, FP time = 0.0, bool status = true) override {
            if (status) return mulRotationMatrix(fields);
            return FP3(0.0, 0.0, 0.0);
        }

        FP3 mulRotationMatrix(const FP3& coords) {
            FP3 directCoords;
            directCoords.x =
                rotationMatrix[0][0] * coords.x +
                rotationMatrix[0][1] * coords.y +
                rotationMatrix[0][2] * coords.z;
            directCoords.y =
                rotationMatrix[1][0] * coords.x +
                rotationMatrix[1][1] * coords.y +
                rotationMatrix[1][2] * coords.z;
            directCoords.z =
                rotationMatrix[2][0] * coords.x +
                rotationMatrix[2][1] * coords.y +
                rotationMatrix[2][2] * coords.z;
            return directCoords;
        }

        FP3 mulInverseRotationMatrix(const FP3& coords) {
            FP3 inverseCoords;
            inverseCoords.x =
                rotationMatrix[0][0] * coords.x +
                rotationMatrix[1][0] * coords.y +
                rotationMatrix[2][0] * coords.z;
            inverseCoords.y =
                rotationMatrix[0][1] * coords.x +
                rotationMatrix[1][1] * coords.y +
                rotationMatrix[2][1] * coords.z;
            inverseCoords.z =
                rotationMatrix[0][2] * coords.x +
                rotationMatrix[1][2] * coords.y +
                rotationMatrix[2][2] * coords.z;
            return inverseCoords;
        }

        Mapping* createInstance() override {
            return new RotationMapping(*this);
        }

        TRotMatrix rotationMatrix = { 0.0 };

    protected:

        void createRotationMatrix(TRotMatrix& matr, int axis, FP angle)
        {
            int axis1 = (axis + 1) % dim, axis2 = (axis + 2) % dim;
            matr[axis1][axis1] = cos(angle);
            matr[axis1][axis2] = -sin(angle);
            matr[axis2][axis1] = sin(angle);
            matr[axis2][axis2] = cos(angle);
            matr[axis][axis] = (FP)1.0;
        }

        void mulMatrices(const TRotMatrix& a, const TRotMatrix& b, TRotMatrix& res) {
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++) {
                    FP sum = 0.0;
                    for (int k = 0; k < dim; k++)
                        sum += a[i][k] * b[k][j];
                    res[i][j] = sum;
                }
        }

    };


    class ShiftMapping : public Mapping {

    public:

        ShiftMapping(const FP3& shift) : shift(shift) {}

        FP3 getDirectCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            setOkStatus(status);
            return coords + shift;
        }

        FP3 getInverseCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            setOkStatus(status);
            return coords - shift;
        }

        Mapping* createInstance() override {
            return new ShiftMapping(*this);
        }

        FP3 shift;

    };


    class ScaleMapping : public Mapping {

    public:

        ScaleMapping(CoordinateEnum axis, FP coef) : axis(axis), coef(coef) {}

        FP3 getDirectCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            setOkStatus(status);
            FP3 directCoords = coords;
            directCoords[(int)axis] *= coef;
            return directCoords;
        }

        FP3 getInverseCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            setOkStatus(status);
            FP3 inverseCoords = coords;
            inverseCoords[(int)axis] /= coef;
            return inverseCoords;
        }

        Mapping* createInstance() override {
            return new ScaleMapping(*this);
        }

        FP coef;
        CoordinateEnum axis;

    };


    class TightFocusingMapping : public Mapping {

    public:

        TightFocusingMapping(FP R0, FP L, FP D, CoordinateEnum axis = CoordinateEnum::x) :
            Rmax(R0 + 0.5*L), ifCut(true),
            periodicalMapping(axis, -R0 - D + 0.5*L, -R0 + 0.5*L)
        {}

        void setIfCut(bool ifCut = true) {
            this->ifCut = ifCut;
        }

        FP3 getDirectCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {
            FP ct = constants::c*time;
            FP r = coords.norm();

            FP3 directCoords = coords;
            if (ifCut) this->setFailStatus(status);
            else setOkStatus(status);

            if (coords[(int)periodicalMapping.axis] < periodicalMapping.cMin + ct ||
                coords[(int)periodicalMapping.axis] >= periodicalMapping.cMax + ct)
                return coords;

            int nPeriods = 0;
            FP shift = 0;

            if (periodicalMapping.cMax + ct < 0) {
                nPeriods = int(-(periodicalMapping.cMin + ct) / periodicalMapping.D) + 1;  // целая часть сверху
                shift = periodicalMapping.D;
            }
            else if (periodicalMapping.cMin + ct > 0) {
                nPeriods = int((periodicalMapping.cMax + ct) / periodicalMapping.D) + 1;  // целая часть сверху
                shift = -periodicalMapping.D;
            }
            else nPeriods = 1;

            FP3 coordsShift = coords;

            for (int i = 0; i < nPeriods; i++) {

                coordsShift[(int)periodicalMapping.axis] += shift;

                if (ifInArea(coordsShift, time)) {
                    directCoords = coordsShift;
                    setOkStatus(status);
                    break;
                }
            }

            return directCoords;
        }

        FP3 getInverseCoords(const FP3& coords, FP time = 0.0, bool* status = 0) override {

            setOkStatus(status);

            if (ifCut && !ifInArea(coords, time)) {
                setFailStatus(status);
            }

            return periodicalMapping.getInverseCoords(coords);

        }

        FP getMinCoord() const { return periodicalMapping.cMin; }
        FP getMaxCoord() const { return periodicalMapping.cMax; }

        bool ifInArea(const FP3& coords, FP time) {
            FP ct = constants::c*time;
            FP r = coords.norm();

            if (periodicalMapping.cMax + ct < 0) {
                if ((r >= Rmax - ct) || (r < -periodicalMapping.cMax - ct) ||
                    (coords[(int)periodicalMapping.axis] > 0))
                    return false;
            }
            else if ((periodicalMapping.cMax + ct >= 0) && (-Rmax + ct <= 0))
            {
                if (coords[(int)periodicalMapping.axis] < 0 && r > periodicalMapping.cMax + Rmax)
                    return false;
                if (coords[(int)periodicalMapping.axis] >= 0 && r > periodicalMapping.cMax + ct)
                    return false;
            }
            else if (-Rmax + ct > 0)
            {
                if ((r <= -Rmax + ct) || (r > periodicalMapping.cMax + ct) ||
                    (coords[(int)periodicalMapping.axis] < 0))
                    return false;
            }

            return true;
        }

        Mapping* createInstance() override {
            return new TightFocusingMapping(*this);
        }

        FP Rmax;
        bool ifCut = true;
        PeriodicalMapping periodicalMapping;

    };

}