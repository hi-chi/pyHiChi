#pragma once

#include "Vectors.h"
#include "VectorsProxy.h"

using namespace std;
using namespace pfc;

namespace pfc {
    struct ValueField
    {
    public:
        ValueField(FP3 E_value, FP3 B_value) :
            E(E_value), B(B_value)
        {}

        ValueField(FP Ex = 0.0, FP Ey = 0.0, FP Ez = 0.0,
            FP Bx = 0.0, FP By = 0.0, FP Bz = 0.0) :
            E(Ex, Ey, Ez), B(Bx, By, Bz)
        {}

        void setE(const FP3 newE)
        { E = newE; }
        FP3 getE()
        { return E; }

        void setB(const FP3 newB)
        { B = newB; }
        FP3 getB()
        { return B; }

    public:
        FP3 E;
        FP3 B;

        friend class ValueFieldProxy;
    };

    /*struct ValueFieldProxy
    {
    public:
        ValueFieldProxy(ValueField& field) :
            E(field.E), B(field.B)
        {}
        ValueFieldProxy(FP3& E_value, FP3& B_value) :
            E(E_value), B(B_value)
        {}

        ValueFieldProxy(FP& Ex, FP& Ey, FP& Ez, FP& Bx, FP& By, FP& Bz) :
            E(Ex, Ey, Ez), B(Bx, By, Bz)
        {}

        void setE(const FP3 newE)
        { E = newE; }
        FP3Proxy getE()
        { return E; }

        void setB(const FP3 newB)
        { B = newB; }
        FP3Proxy getB()
        { return B; }

    private:
        FP3Proxy E;
        FP3Proxy B;
    };*/

}