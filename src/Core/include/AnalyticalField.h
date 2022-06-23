#pragma once
#include "Vectors.h"


namespace pfc {

    class AnalyticalField {

    public:

        typedef FP(*TFunc)(FP, FP, FP, FP);

        AnalyticalField(FP dt) : dt(dt), globalTime(0.0) {}

        void setTimeStep(FP dt) { this->dt = dt; }
        void updateFields() { this->globalTime += dt; }

        void setTime(FP t) { this->globalTime = t; }
        FP getTime() { return this->globalTime; }

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

        FP3 getE(const FP3& coord) const {
            return this->getE(coord.x, coord.y, coord.z, this->globalTime);
        }
        FP3 getB(const FP3& coord) const {
            return this->getB(coord.x, coord.y, coord.z, this->globalTime);
        }
        FP3 getJ(const FP3& coord) const {
            return this->getJ(coord.x, coord.y, coord.z, this->globalTime);
        }

        FP getEx(const FP3& coord) const { return this->funcEx(coord.x, coord.y, coord.z, this->globalTime); }
        FP getEy(const FP3& coord) const { return this->funcEy(coord.x, coord.y, coord.z, this->globalTime); }
        FP getEz(const FP3& coord) const { return this->funcEz(coord.x, coord.y, coord.z, this->globalTime); }

        FP getBx(const FP3& coord) const { return this->funcBx(coord.x, coord.y, coord.z, this->globalTime); }
        FP getBy(const FP3& coord) const { return this->funcBy(coord.x, coord.y, coord.z, this->globalTime); }
        FP getBz(const FP3& coord) const { return this->funcBz(coord.x, coord.y, coord.z, this->globalTime); }

        FP getJx(const FP3& coord) const { return this->funcJx(coord.x, coord.y, coord.z, this->globalTime); }
        FP getJy(const FP3& coord) const { return this->funcJy(coord.x, coord.y, coord.z, this->globalTime); }
        FP getJz(const FP3& coord) const { return this->funcJz(coord.x, coord.y, coord.z, this->globalTime); }

        FP globalTime;
        FP dt;
        FP timeShiftE = 0.0, timeShiftB = 0.0, timeShiftJ = 0.0;

    private:

        TFunc funcEx, funcEy, funcEz, funcBx, funcBy, funcBz, funcJx, funcJy, funcJz;

    };

}