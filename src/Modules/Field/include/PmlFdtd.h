#pragma once
#include <limits>

#include "Pml.h"

namespace pfc {

    class PmlFdtd : public PmlReal<YeeGrid>
    {
    public:
        PmlFdtd(YeeGrid* grid, FP dt, Int3 sizePML) :
            PmlReal(grid, dt, sizePML) {}

        void updateB();
        void updateE();

    private:

        void updateB3D();
        void updateB2D();
        void updateB1D();
        void updateE3D();
        void updateE2D();
        void updateE1D();
    };

    inline void PmlFdtd::updateB()
    {
        if (this->grid->dimensionality == 3)
            updateB3D();
        else if (this->grid->dimensionality == 2)
            updateB2D();
        else if (this->grid->dimensionality == 1)
            updateB1D();
    }

    inline void PmlFdtd::updateB3D()
    {
        const int numPmlIndicesB = this->bIndex.size();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < numPmlIndicesB; ++idx)
        {
            int i = bIndex[idx].x, j = bIndex[idx].y, k = bIndex[idx].z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: By(i, ..., ...), Bz(i, ..., ...) -> sigmaX(i)
            sigma = this->computeSigma(this->grid->ByPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = cdt / dx;

            byx[idx] = coeff1 * byx[idx] + coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i - 1, j, k));
            bzx[idx] = coeff1 * bzx[idx] - coeff2 *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i - 1, j, k));

            // Yee Grid indices: Bx(..., j, ...), Bz(..., j, ...) -> sigmaY(j)
            sigma = this->computeSigma(this->grid->BzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dy);
            else coeff2 = cdt / dy;

            bxy[idx] = coeff1 * bxy[idx] - coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i, j - 1, k));
            bzy[idx] = coeff1 * bzy[idx] + coeff2 *
                (this->grid->Ex(i, j, k) - this->grid->Ex(i, j - 1, k));

            // Yee Grid indices: Bx(..., ..., k), By(..., ..., k) -> sigmaZ(k)
            sigma = this->computeSigma(this->grid->BxPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dz);
            else coeff2 = cdt / dz;

            bxz[idx] = coeff1 * bxz[idx] + coeff2 *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i, j, k - 1));
            byz[idx] = coeff1 * byz[idx] - coeff2 *
                (this->grid->Ex(i, j, k) - this->grid->Ex(i, j, k - 1));

            this->grid->Bx(i, j, k) = bxy[idx] + bxz[idx];
            this->grid->By(i, j, k) = byx[idx] + byz[idx];
            this->grid->Bz(i, j, k) = bzx[idx] + bzy[idx];
        }
    }

    inline void PmlFdtd::updateB2D()
    {
        const int numPmlIndicesB = this->bIndex.size();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < numPmlIndicesB; ++idx)
        {
            int i = bIndex[idx].x, j = bIndex[idx].y, k = bIndex[idx].z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: By(i, ..., ...), Bz(i, ..., ...) -> sigmaX(i)
            sigma = this->computeSigma(this->grid->ByPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = cdt / dx;

            byx[idx] = coeff1 * byx[idx] + coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i - 1, j, k));
            bzx[idx] = coeff1 * bzx[idx] - coeff2 *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i - 1, j, k));

            // Yee Grid indices: Bx(..., j, ...), Bz(..., j, ...) -> sigmaY(j)
            sigma = this->computeSigma(this->grid->BzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dy);
            else coeff2 = cdt / dy;

            bxy[idx] = coeff1 * bxy[idx] - coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i, j - 1, k));
            bzy[idx] = coeff1 * bzy[idx] + coeff2 *
                (this->grid->Ex(i, j, k) - this->grid->Ex(i, j - 1, k));

            // Yee Grid indices: Bx(..., ..., k), By(..., ..., k) -> sigmaZ(k)
            sigma = this->computeSigma(this->grid->BxPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);

            bxz[idx] = coeff1 * bxz[idx];
            byz[idx] = coeff1 * byz[idx];

            this->grid->Bx(i, j, k) = bxy[idx] + bxz[idx];
            this->grid->By(i, j, k) = byx[idx] + byz[idx];
            this->grid->Bz(i, j, k) = bzx[idx] + bzy[idx];
        }
    }

    inline void PmlFdtd::updateB1D()
    {
        const int numPmlIndicesB = this->bIndex.size();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < numPmlIndicesB; ++idx)
        {
            int i = bIndex[idx].x, j = bIndex[idx].y, k = bIndex[idx].z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: By(i, ..., ...), Bz(i, ..., ...) -> sigmaX(i)
            sigma = this->computeSigma(this->grid->ByPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = cdt / dx;

            byx[idx] = coeff1 * byx[idx] + coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i - 1, j, k));
            bzx[idx] = coeff1 * bzx[idx] - coeff2 *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i - 1, j, k));

            // Yee Grid indices: Bx(..., j, ...), Bz(..., j, ...) -> sigmaY(j)
            sigma = this->computeSigma(this->grid->BzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);

            bxy[idx] = coeff1 * bxy[idx];
            bzy[idx] = coeff1 * bzy[idx];

            // Yee Grid indices: Bx(..., ..., k), By(..., ..., k) -> sigmaZ(k)
            sigma = this->computeSigma(this->grid->BxPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);

            bxz[idx] = coeff1 * bxz[idx];
            byz[idx] = coeff1 * byz[idx];

            this->grid->Bx(i, j, k) = bxy[idx] + bxz[idx];
            this->grid->By(i, j, k) = byx[idx] + byz[idx];
            this->grid->Bz(i, j, k) = bzx[idx] + bzy[idx];
        }
    }

    inline void PmlFdtd::updateE()
    {
        if (this->grid->dimensionality == 3)
            updateE3D();
        else if (this->grid->dimensionality == 2)
            updateE2D();
        else if (this->grid->dimensionality == 1)
            updateE1D();
    }

    inline void PmlFdtd::updateE3D()
    {
        const int numPmlIndicesE = this->eIndex.size();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < numPmlIndicesE; ++idx)
        {
            int i = eIndex[idx].x, j = eIndex[idx].y, k = eIndex[idx].z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: Ey(i+1/2, ..., ...), Ez(i+1/2, ..., ...) -> sigmaX(i+1/2)
            sigma = this->computeSigma(this->grid->EyPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = -cdt / dx;

            eyx[idx] = coeff1 * eyx[idx] + coeff2 *
                (this->grid->Bz(i + 1, j, k) - this->grid->Bz(i, j, k));
            ezx[idx] = coeff1 * ezx[idx] - coeff2 *
                (this->grid->By(i + 1, j, k) - this->grid->By(i, j, k));

            // Yee Grid indices: Ex(..., j+1/2, ...), Ez(..., j+1/2, ...)-> sigmaY(j+1/2)
            sigma = this->computeSigma(this->grid->EzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dy);
            else coeff2 = -cdt / dy;

            exy[idx] = coeff1 * exy[idx] - coeff2 *
                (this->grid->Bz(i, j + 1, k) - this->grid->Bz(i, j, k));
            ezy[idx] = coeff1 * ezy[idx] + coeff2 *
                (this->grid->Bx(i, j + 1, k) - this->grid->Bx(i, j, k));

            // Yee Grid indices: Ex(..., ..., k+1/2), Ey(..., ..., k+1/2)-> sigmaZ(k+1/2)
            sigma = this->computeSigma(this->grid->ExPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dz);
            else coeff2 = -cdt / dz;

            exz[idx] = coeff1 * exz[idx] + coeff2 *
                (this->grid->By(i, j, k + 1) - this->grid->By(i, j, k));
            eyz[idx] = coeff1 * eyz[idx] - coeff2 *
                (this->grid->Bx(i, j, k + 1) - this->grid->Bx(i, j, k));

            this->grid->Ex(i, j, k) = exy[idx] + exz[idx];
            this->grid->Ey(i, j, k) = eyx[idx] + eyz[idx];
            this->grid->Ez(i, j, k) = ezx[idx] + ezy[idx];
        }
    }

    inline void PmlFdtd::updateE2D()
    {
        const int numPmlIndicesE = this->eIndex.size();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < numPmlIndicesE; ++idx)
        {
            int i = eIndex[idx].x, j = eIndex[idx].y, k = eIndex[idx].z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: Ey(i+1/2, ..., ...), Ez(i+1/2, ..., ...) -> sigmaX(i+1/2)
            sigma = this->computeSigma(this->grid->EyPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = -cdt / dx;

            eyx[idx] = coeff1 * eyx[idx] + coeff2 *
                (this->grid->Bz(i + 1, j, k) - this->grid->Bz(i, j, k));
            ezx[idx] = coeff1 * ezx[idx] - coeff2 *
                (this->grid->By(i + 1, j, k) - this->grid->By(i, j, k));

            // Yee Grid indices: Ex(..., j+1/2, ...), Ez(..., j+1/2, ...)-> sigmaY(j+1/2)
            sigma = this->computeSigma(this->grid->EzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dy);
            else coeff2 = -cdt / dy;

            exy[idx] = coeff1 * exy[idx] - coeff2 *
                (this->grid->Bz(i, j + 1, k) - this->grid->Bz(i, j, k));
            ezy[idx] = coeff1 * ezy[idx] + coeff2 *
                (this->grid->Bx(i, j + 1, k) - this->grid->Bx(i, j, k));

            // Yee Grid indices: Ex(..., ..., k+1/2), Ey(..., ..., k+1/2)-> sigmaZ(k+1/2)
            sigma = this->computeSigma(this->grid->ExPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);

            exz[idx] = coeff1 * exz[idx];
            eyz[idx] = coeff1 * eyz[idx];

            this->grid->Ex(i, j, k) = exy[idx] + exz[idx];
            this->grid->Ey(i, j, k) = eyx[idx] + eyz[idx];
            this->grid->Ez(i, j, k) = ezx[idx] + ezy[idx];
        }
    }

    inline void PmlFdtd::updateE1D()
    {
        const int numPmlIndicesE = this->eIndex.size();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < numPmlIndicesE; ++idx)
        {
            int i = eIndex[idx].x, j = eIndex[idx].y, k = eIndex[idx].z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: Ey(i+1/2, ..., ...), Ez(i+1/2, ..., ...) -> sigmaX(i+1/2)
            sigma = this->computeSigma(this->grid->EyPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = -cdt / dx;

            eyx[idx] = coeff1 * eyx[idx] + coeff2 *
                (this->grid->Bz(i + 1, j, k) - this->grid->Bz(i, j, k));
            ezx[idx] = coeff1 * ezx[idx] + coeff2 *
                (this->grid->By(i + 1, j, k) - this->grid->By(i, j, k));

            // Yee Grid indices: Ex(..., j+1/2, ...), Ez(..., j+1/2, ...)-> sigmaY(j+1/2)
            sigma = this->computeSigma(this->grid->EzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);

            exy[idx] = coeff1 * exy[idx];
            ezy[idx] = coeff1 * ezy[idx];

            // Yee Grid indices: Ex(..., ..., k+1/2), Ey(..., ..., k+1/2)-> sigmaZ(k+1/2)
            sigma = this->computeSigma(this->grid->ExPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);

            exz[idx] = coeff1 * exz[idx];
            eyz[idx] = coeff1 * eyz[idx];

            this->grid->Ex(i, j, k) = exy[idx] + exz[idx];
            this->grid->Ey(i, j, k) = eyx[idx] + eyz[idx];
            this->grid->Ez(i, j, k) = ezx[idx] + ezy[idx];
        }
    }
}
