#pragma once
#include "Vectors.h"

namespace pfc {

    class Mapping {

    public:

        // create mapping on a segment [a,b) for grid
        Mapping() {}

        // return direct coords, write to status true if it is possible
        virtual FP3 getDirectCoords(const FP3& coords, bool* status = 0) {
            *status = true;
            return coords;
        }

        // return inverse coords, write to status true if it is possible
        virtual FP3 getInverseCoords(const FP3& coords, bool* status = 0) {
            *status = true;
            return coords;
        }

    };


    class IdentityMapping : public Mapping {

    public:

        // create identity mapping on a segment [a,b)
        IdentityMapping(const FP3& a, const FP3& b): a(a), b(b) {}

        FP3 getDirectCoords(const FP3& coords, bool* status = 0) override {
            *status = (coords >= a && coords < b) ? true : false;
            return coords;
        }

        FP3 getInverseCoords(const FP3& coords, bool* status = 0) override {
            *status = (coords >= a && coords < b) ? true : false;
            return coords;
        }

        FP3 a, b;

    };

}