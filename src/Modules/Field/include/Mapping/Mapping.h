#pragma once
#include "Vectors.h"
#include <functional>

namespace pfc {

    class Mapping {

    public:

        Mapping() : time(0.0) {}

		virtual FP3 getDirectCoords(const FP3& coords, bool* status = 0) {
			setFailStatus(status);
			return coords;
		}

		virtual FP3 getInverseCoords(const FP3& coords, bool* status = 0) {
			setFailStatus(status);
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

        FP time = 0.0;

    };


    class IdentityMapping : public Mapping {

    public:

        // create identity mapping on a segment [a,b)
		IdentityMapping(FP3 a, FP3 b) : a(a), b(b) {}

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
		PeriodicalXMapping(FP xMin, FP xMax) : xMin(xMin), xMax(xMax), D(xMax-xMin) {}

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

}