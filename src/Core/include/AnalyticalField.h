#pragma once
#include "Vectors.h"


namespace pfc {

    class AnalyticalField {

    public:

        typedef FP(*TFunc)(FP, FP, FP, FP);

        AnalyticalField(FP dt) : dt(dt), globalTime(0.0) {}


        void setTimeStep(double dt) { this->dt = dt; }
        void updateFields() { this->globalTime += dt; }


        void setE(TFunc funcEx, TFunc funcEy, TFunc funcEz) { 
            this->funcEx = funcEx;
            this->funcEy = funcEy;
            this->funcEz = funcEz;
        }
        
        void setB(TFunc funcBx, TFunc funcBy, TFunc funcBz) {
            this->funcBx = funcBx;
            this->funcBy = funcBy;
            this->funcBz = funcBz;
        }

        void setJ(TFunc funcJx, TFunc funcJy, TFunc funcJz) {
            this->funcJx = funcJx;
            this->funcJy = funcJy;
            this->funcJz = funcJz;
        }


        FP3 getE(FP x, FP y, FP z, FP t) const {
            return FP3(this->funcEx(x, y, z, t),
                this->funcEy(x, y, z, t),
                this->funcEz(x, y, z, t));
        }

        FP3 getB(FP x, FP y, FP z, FP t) const {
            return FP3(this->funcBx(x, y, z, t),
                this->funcBy(x, y, z, t),
                this->funcBz(x, y, z, t));
        }

        FP3 getJ(FP x, FP y, FP z, FP t) const {
            return FP3(this->funcJx(x, y, z, t),
                this->funcJy(x, y, z, t),
                this->funcJz(x, y, z, t));
        }


        FP3 getE(FP3 coord) const {
            return FP3(this->funcEx(coord.x, coord.y, coord.z, this->globalTime),
                this->funcEy(coord.x, coord.y, coord.z, this->globalTime),
                this->funcEz(coord.x, coord.y, coord.z, this->globalTime));
        }

        FP3 getB(FP3 coord) const {
            return FP3(this->funcBx(coord.x, coord.y, coord.z, this->globalTime),
                this->funcBy(coord.x, coord.y, coord.z, this->globalTime),
                this->funcBz(coord.x, coord.y, coord.z, this->globalTime));
        }

        FP3 getJ(FP3 coord) const {
            return FP3(this->funcJx(coord.x, coord.y, coord.z, this->globalTime),
                this->funcJy(coord.x, coord.y, coord.z, this->globalTime),
                this->funcJz(coord.x, coord.y, coord.z, this->globalTime));
        }


        FP globalTime;
        FP dt;

    private:

        TFunc funcEx, funcEy, funcEz, funcBx, funcBy, funcBz, funcJx, funcJy, funcJz;

    };

}