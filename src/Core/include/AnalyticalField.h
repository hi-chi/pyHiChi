#pragma once
#include "Vectors.h"

#include <functional>


namespace pfc {

    namespace analytical_field {
        FP defaultFieldFunction(FP x, FP y, FP z, FP t) {
            return (FP)0.0;
        }
    };

    // Implements general grid methods for analytical field
    class AnalyticalField
    {
    public:

        using FunctionType = std::function<FP(FP, FP, FP, FP)>;

        AnalyticalField()
        {
            FunctionType dfunc = analytical_field::defaultFieldFunction;
            setE(dfunc, dfunc, dfunc);
            setB(dfunc, dfunc, dfunc);
            setJ(dfunc, dfunc, dfunc);
        }

        AnalyticalField(FunctionType funcEx, FunctionType funcEy, FunctionType funcEz,
            FunctionType funcBx, FunctionType funcBy, FunctionType funcBz,
            FunctionType funcJx, FunctionType funcJy, FunctionType funcJz)
        {
            setE(funcEx, funcEy, funcEz);
            setB(funcBx, funcBy, funcBz);
            setJ(funcJx, funcJy, funcJz);
        }

        AnalyticalField(FunctionType funcEx, FunctionType funcEy, FunctionType funcEz,
            FunctionType funcBx, FunctionType funcBy, FunctionType funcBz)
        {
            setE(funcEx, funcEy, funcEz);
            setB(funcBx, funcBy, funcBz);
            FunctionType dfunc = analytical_field::defaultFieldFunction;
            setJ(dfunc, dfunc, dfunc);
        }

        void setE(FunctionType funcEx, FunctionType funcEy, FunctionType funcEz) { 
            this->funcEx = funcEx;
            this->funcEy = funcEy;
            this->funcEz = funcEz;
        }      
        void setB(FunctionType funcBx, FunctionType funcBy, FunctionType funcBz) {
            this->funcBx = funcBx;
            this->funcBy = funcBy;
            this->funcBz = funcBz;
        }
        void setJ(FunctionType funcJx, FunctionType funcJy, FunctionType funcJz) {
            this->funcJx = funcJx;
            this->funcJy = funcJy;
            this->funcJz = funcJz;
        }

        FP3 getE(FP x, FP y, FP z, FP t) const {
            return FP3(
                this->funcEx(x, y, z, t),
                this->funcEy(x, y, z, t),
                this->funcEz(x, y, z, t));
        }
        FP3 getB(FP x, FP y, FP z, FP t) const {
            return FP3(
                this->funcBx(x, y, z, t),
                this->funcBy(x, y, z, t),
                this->funcBz(x, y, z, t));
        }
        FP3 getJ(FP x, FP y, FP z, FP t) const {
            return FP3(
                this->funcJx(x, y, z, t),
                this->funcJy(x, y, z, t),
                this->funcJz(x, y, z, t));
        }

        FP getEx(FP x, FP y, FP z, FP t) const { return this->funcEx(x, y, z, t); }
        FP getEy(FP x, FP y, FP z, FP t) const { return this->funcEy(x, y, z, t); }
        FP getEz(FP x, FP y, FP z, FP t) const { return this->funcEz(x, y, z, t); }        
        FP getBx(FP x, FP y, FP z, FP t) const { return this->funcBx(x, y, z, t); }
        FP getBy(FP x, FP y, FP z, FP t) const { return this->funcBy(x, y, z, t); }
        FP getBz(FP x, FP y, FP z, FP t) const { return this->funcBz(x, y, z, t); }        
        FP getJx(FP x, FP y, FP z, FP t) const { return this->funcJx(x, y, z, t); }
        FP getJy(FP x, FP y, FP z, FP t) const { return this->funcJy(x, y, z, t); }
        FP getJz(FP x, FP y, FP z, FP t) const { return this->funcJz(x, y, z, t); }

        FP3 getE(const FP3& coords) const {
            return FP3(
                this->funcEx(coords.x, coords.y, coords.z, this->globalTime),
                this->funcEy(coords.x, coords.y, coords.z, this->globalTime),
                this->funcEz(coords.x, coords.y, coords.z, this->globalTime));
        }
        FP3 getB(const FP3& coords) const {
            return FP3(
                this->funcBx(coords.x, coords.y, coords.z, this->globalTime),
                this->funcBy(coords.x, coords.y, coords.z, this->globalTime),
                this->funcBz(coords.x, coords.y, coords.z, this->globalTime));
        }
        FP3 getJ(const FP3& coords) const {
            return FP3(
                this->funcJx(coords.x, coords.y, coords.z, this->globalTime),
                this->funcJy(coords.x, coords.y, coords.z, this->globalTime),
                this->funcJz(coords.x, coords.y, coords.z, this->globalTime));
        }

        FP getEx(const FP3& coords) const { return this->funcEx(coords.x, coords.y, coords.z, this->globalTime); }
        FP getEy(const FP3& coords) const { return this->funcEy(coords.x, coords.y, coords.z, this->globalTime); }
        FP getEz(const FP3& coords) const { return this->funcEz(coords.x, coords.y, coords.z, this->globalTime); }
        FP getBx(const FP3& coords) const { return this->funcBx(coords.x, coords.y, coords.z, this->globalTime); }
        FP getBy(const FP3& coords) const { return this->funcBy(coords.x, coords.y, coords.z, this->globalTime); }
        FP getBz(const FP3& coords) const { return this->funcBz(coords.x, coords.y, coords.z, this->globalTime); }
        FP getJx(const FP3& coords) const { return this->funcJx(coords.x, coords.y, coords.z, this->globalTime); }
        FP getJy(const FP3& coords) const { return this->funcJy(coords.x, coords.y, coords.z, this->globalTime); }
        FP getJz(const FP3& coords) const { return this->funcJz(coords.x, coords.y, coords.z, this->globalTime); }

        void save(std::ostream& ostr) {
            ostr.write((char*)&globalTime, sizeof(globalTime));
            // TODO: save functions
        }
        void load(std::istream& istr) {
            istr.read((char*)&globalTime, sizeof(globalTime));
            // TODO: load functions
        }

        FP globalTime = 0.0;

    private:

        FunctionType funcEx, funcEy, funcEz, funcBx, funcBy, funcBz, funcJx, funcJy, funcJz;

    };

}