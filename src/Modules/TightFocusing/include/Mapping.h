#pragma once
#include "Vectors.h"
#include <functional>

namespace pfc {

    class Mapping {

    public:

        Mapping() : time(0) {}

        // return direct coords, write to status true if it is possible
        virtual FP3 getDirectCoords(const FP3& coords, bool* status = 0) {
            if (status) *status = true;
            return coords;
        }

        // return inverse coords, write to status true if it is possible
        virtual FP3 getInverseCoords(const FP3& coords, bool* status = 0) {
            if (status) *status = true;
            return coords;
        }

        void advanceTime(FP timeStep) {
            time += timeStep;
        }

		void setTime(FP time) {
			this->time = time;
		}

    protected:

        void setFailStatus(bool* status) {
            if (status) *status = false;
        }
        void setOkStatus(bool* status) {
            if (status) *status = true;
        }

        FP time = 0;

    };


    class IdentityMapping : public Mapping {

    public:

        // create identity mapping on a segment [a,b)
        IdentityMapping(const FP3& a, const FP3& b): a(a), b(b) {}

        FP3 getDirectCoords(const FP3& coords, bool* status = 0) override {
            if (status) *status = (coords >= a && coords < b) ? true : false;
            return coords;
        }

        FP3 getInverseCoords(const FP3& coords, bool* status = 0) override {
            if (status) *status = (coords >= a && coords < b) ? true : false;
            return coords;
        }

        FP3 a, b;

    };


    class PeriodicalXMapping : public Mapping {

    public:

        // create periodical mapping: x = ...[xMin, xMax)[xMin, xMax)[xMin, xMax)...
        PeriodicalXMapping(FP xMin, FP xMax) : xMin(xMin), xMax(xMax), D(xMax-xMin){}

        FP3 getDirectCoords(const FP3& coords, bool* status = 0) override {
            if (status) *status = (coords.x >= xMin && coords.x < xMax) ? true : false;
            return coords;
        }

        FP3 getInverseCoords(const FP3& coords, bool* status = 0) override {
            setOkStatus(status);
            FP3 inverseCoords = coords;
            double tmp;
            FP frac = std::modf((coords.x - xMin) / D, &tmp);
            inverseCoords.x = xMin + (frac >= 0 ? frac : (1 + frac)) * D;
            return inverseCoords;
        }

        FP xMin, xMax, D;

    };

    class TightFocusingMapping : public PeriodicalXMapping {

    public:

        TightFocusingMapping(FP R0, FP L, FP D, FP cutAngle = 0.5*constants::pi) :
            PeriodicalXMapping(-R0 - D + 0.5*L, -R0 + 0.5*L),
            xL(-R0 - 0.5*L), cutAngle(cutAngle) {}

        FP3 getDirectCoords(const FP3& coords, bool* status = 0) override {
            FP ct = constants::c*time;
            FP r = coords.norm();

            FP3 directCoords = coords;
            setFailStatus(status);

            if (coords.x < xMin + ct || coords.x >= xMax + ct)
                return coords;

            int nPeriods = 0;
            int shiftSign = 0;

            if (xMax + ct < 0) {
                nPeriods = int(((xMax + ct)*cos(cutAngle) - (xMin + ct)) / D) + 1;  // целая часть сверху
                shiftSign = 1;
            }
            else if (xMin + ct > 0) {
                nPeriods = int(((xMax + ct) - (xMin + ct)*cos(cutAngle)) / D) + 1;  // целая часть сверху
                shiftSign = -1;
            }
            else nPeriods = 1;

            for (int i = 0; i < nPeriods; i++) {
                FP shift = i * D * shiftSign;

                FP3 coordsShift(coords.x + shift, coords.y, coords.z);
                FP rShift = coordsShift.norm();

                if (ifInArea(coordsShift)) {
                    directCoords = coordsShift;
                    setOkStatus(status);
                    break;
                }
            }

            return directCoords;
        }

        FP3 getInverseCoords(const FP3& coords, bool* status = 0) override {
            
            if (ifInArea(coords)) {
                setOkStatus(status);
            }
            else {
                setFailStatus(status);
            }

            return PeriodicalXMapping::getInverseCoords(coords);

        }

        FP getxMin() const { return xMin; }
        FP getxMax() const { return xMax; }

        bool ifInArea(const FP3& coords) {
            FP ct = constants::c*time;
            FP r = coords.norm();
            FP angle = atan(abs(sqrt(coords.y*coords.y + coords.z*coords.z)/coords.x));

            if (xMax + ct < 0) {
                if ((r >= -xL - ct) || (r < -xMax - ct) || (coords.x > 0))
                    return false;
                if (angle > cutAngle)
                    return false;
            }
            else if ((xMax + ct >= 0) && (xL + ct <= 0))
            {
                if (coords.x < 0 && r > xMax - xL)
                    return false;
                if (coords.x >= 0 && r > xMax + ct)
                    return false;
            }
            else if (xL + ct > 0)
            {
                if ((r <= xL + ct) || (r > xMax + ct) || (coords.x < 0))
                    return false;
                if (angle > cutAngle)
                    return false;
            }

            return true;
        }

        FP cutAngle, xL;

    };

}