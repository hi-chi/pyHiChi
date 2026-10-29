#pragma once
#include <cmath>

#include "FP.h"
#include "Vectors.h"

namespace pfc {
    inline FP formfactorTSC(FP x)
    {
        FP xa = fabs(x);
#ifdef __MIC
        return (xa <= (FP)0.5) * ((FP)0.75 - x * x) + (x > (FP)0.5) * ((FP)0.5 * ((FP)1.5 - xa) * ((FP)1.5 - xa));
#else
        if (xa <= (FP)0.5)
            return (FP)0.75 - x * x;
        else
            return (FP)0.5 * ((FP)1.5 - xa) * ((FP)1.5 - xa);
#endif
    }

    inline FP formfactorPCS(FP x)
    {
        FP xa = fabs(x);
#ifdef __MIC
        return (xa <= (FP)1.0) * (((FP)4.0 - 6 * xa * xa + 3 * xa * xa * xa) / (FP)6.0)
            + (xa > (FP)1.0) * (((FP)2.0 - xa) * ((FP)2.0 - xa) * ((FP)2.0 - xa) / (FP)6.0);
#else
        if (xa <= (FP)1.0)
            return ((FP)4.0 - 6 * xa * xa + 3 * xa * xa * xa) / (FP)6.0;
        else
            return ((FP)2.0 - xa) * ((FP)2.0 - xa) * ((FP)2.0 - xa) / (FP)6.0;
#endif
    }

    inline void formfactorFourthOrder(FP coeff, FP c[])
    {
        c[0] = (coeff + (FP)1) * coeff * (coeff - (FP)1) * (coeff - (FP)2) / (FP)24;
        c[1] = -(coeff + (FP)2) * coeff * (coeff - (FP)1) * (coeff - (FP)2) / (FP)6;
        c[2] = (coeff + (FP)2) * (coeff + (FP)1) * (coeff - (FP)1) * (coeff - (FP)2) / (FP)4;
        c[3] = -(coeff + (FP)2) * (coeff + (FP)1) * coeff * (coeff - (FP)2) / (FP)6;
        c[4] = (coeff + (FP)2) * (coeff + (FP)1) * coeff * (coeff - (FP)1) / (FP)24;
    }

    class FormFactor {
    public:
        void operator()(const FP3& coords) {};
    };

    class FormFactorCIC : public FormFactor {
    public:
        void operator()(FP3 coords) {
            for (int i = 0; i < 3; ++i) {
                c[i][0] = (FP)1 - coords[i];
                c[i][1] = coords[i];
            }
        }

        FP c[3][2];
    };

    class FormFactorTSC : public FormFactor {
    public:
        void operator()(const FP3& coords) {
            for (int ii = 0; ii < 3; ii++)
                c[0][ii] = formfactorTSC(FP(ii - 1) - coords.x);
            for (int jj = 0; jj < 3; jj++)
                c[1][jj] = formfactorTSC(FP(jj - 1) - coords.y);
            for (int kk = 0; kk < 3; kk++)
                c[2][kk] = formfactorTSC(FP(kk - 1) - coords.z);
        }

        FP c[3][3];
    };
}
