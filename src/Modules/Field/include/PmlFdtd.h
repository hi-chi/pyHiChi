#pragma once
#include <limits>

#include "Pml.h"

namespace pfc {

    class PmlFdtd : public PmlReal<YeeGrid>
    {
    public:

        PmlFdtd(YeeGrid* grid, FP dt, Int3 domainIndexBegin, Int3 domainIndexEnd,
            Int3 sizePML, FP nPmlParam = (FP)4.0, FP r0PmlParam = (FP)1e-8) :
            PmlReal(grid, dt, domainIndexBegin, domainIndexEnd, sizePML, nPmlParam, r0PmlParam)
        {}

        // constructor for loading
        explicit PmlFdtd(YeeGrid* grid, FP dt, Int3 domainIndexBegin, Int3 domainIndexEnd) :
            PmlReal(grid, dt, domainIndexBegin, domainIndexEnd)
        {}

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
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: By(i, ..., ...), Bz(i, ..., ...) -> sigmaX(i)
            sigma = this->computeSigma(this->grid->ByPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = cdt / dx;

            this->splitGrid->byx[idx] = coeff1 * this->splitGrid->byx[idx] + coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i - 1, j, k));
            this->splitGrid->bzx[idx] = coeff1 * this->splitGrid->bzx[idx] - coeff2 *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i - 1, j, k));

            // Yee Grid indices: Bx(..., j, ...), Bz(..., j, ...) -> sigmaY(j)
            sigma = this->computeSigma(this->grid->BzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dy);
            else coeff2 = cdt / dy;

            this->splitGrid->bxy[idx] = coeff1 * this->splitGrid->bxy[idx] - coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i, j - 1, k));
            this->splitGrid->bzy[idx] = coeff1 * this->splitGrid->bzy[idx] + coeff2 *
                (this->grid->Ex(i, j, k) - this->grid->Ex(i, j - 1, k));

            // Yee Grid indices: Bx(..., ..., k), By(..., ..., k) -> sigmaZ(k)
            sigma = this->computeSigma(this->grid->BxPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dz);
            else coeff2 = cdt / dz;

            this->splitGrid->bxz[idx] = coeff1 * this->splitGrid->bxz[idx] + coeff2 *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i, j, k - 1));
            this->splitGrid->byz[idx] = coeff1 * this->splitGrid->byz[idx] - coeff2 *
                (this->grid->Ex(i, j, k) - this->grid->Ex(i, j, k - 1));

            this->grid->Bx(i, j, k) = this->splitGrid->bxy[idx] + this->splitGrid->bxz[idx];
            this->grid->By(i, j, k) = this->splitGrid->byx[idx] + this->splitGrid->byz[idx];
            this->grid->Bz(i, j, k) = this->splitGrid->bzx[idx] + this->splitGrid->bzy[idx];
        }
    }

    inline void PmlFdtd::updateB2D()
    {
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: By(i, ..., ...), Bz(i, ..., ...) -> sigmaX(i)
            sigma = this->computeSigma(this->grid->ByPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = cdt / dx;

            this->splitGrid->byx[idx] = coeff1 * this->splitGrid->byx[idx] + coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i - 1, j, k));
            this->splitGrid->bzx[idx] = coeff1 * this->splitGrid->bzx[idx] - coeff2 *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i - 1, j, k));

            // Yee Grid indices: Bx(..., j, ...), Bz(..., j, ...) -> sigmaY(j)
            sigma = this->computeSigma(this->grid->BzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dy);
            else coeff2 = cdt / dy;

            this->splitGrid->bxy[idx] = coeff1 * this->splitGrid->bxy[idx] - coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i, j - 1, k));
            this->splitGrid->bzy[idx] = coeff1 * this->splitGrid->bzy[idx] + coeff2 *
                (this->grid->Ex(i, j, k) - this->grid->Ex(i, j - 1, k));

            // Yee Grid indices: Bx(..., ..., k), By(..., ..., k) -> sigmaZ(k)
            sigma = this->computeSigma(this->grid->BxPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);

            this->splitGrid->bxz[idx] = coeff1 * this->splitGrid->bxz[idx];
            this->splitGrid->byz[idx] = coeff1 * this->splitGrid->byz[idx];

            this->grid->Bx(i, j, k) = this->splitGrid->bxy[idx] + this->splitGrid->bxz[idx];
            this->grid->By(i, j, k) = this->splitGrid->byx[idx] + this->splitGrid->byz[idx];
            this->grid->Bz(i, j, k) = this->splitGrid->bzx[idx] + this->splitGrid->bzy[idx];
        }
    }

    inline void PmlFdtd::updateB1D()
    {
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: By(i, ..., ...), Bz(i, ..., ...) -> sigmaX(i)
            sigma = this->computeSigma(this->grid->ByPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = -(coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = cdt / dx;

            this->splitGrid->byx[idx] = coeff1 * this->splitGrid->byx[idx] + coeff2 *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i - 1, j, k));
            this->splitGrid->bzx[idx] = coeff1 * this->splitGrid->bzx[idx] - coeff2 *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i - 1, j, k));

            // Yee Grid indices: Bx(..., j, ...), Bz(..., j, ...) -> sigmaY(j)
            sigma = this->computeSigma(this->grid->BzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);

            this->splitGrid->bxy[idx] = coeff1 * this->splitGrid->bxy[idx];
            this->splitGrid->bzy[idx] = coeff1 * this->splitGrid->bzy[idx];

            // Yee Grid indices: Bx(..., ..., k), By(..., ..., k) -> sigmaZ(k)
            sigma = this->computeSigma(this->grid->BxPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);

            this->splitGrid->bxz[idx] = coeff1 * this->splitGrid->bxz[idx];
            this->splitGrid->byz[idx] = coeff1 * this->splitGrid->byz[idx];

            this->grid->Bx(i, j, k) = this->splitGrid->bxy[idx] + this->splitGrid->bxz[idx];
            this->grid->By(i, j, k) = this->splitGrid->byx[idx] + this->splitGrid->byz[idx];
            this->grid->Bz(i, j, k) = this->splitGrid->bzx[idx] + this->splitGrid->bzy[idx];
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
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: Ey(i+1/2, ..., ...), Ez(i+1/2, ..., ...) -> sigmaX(i+1/2)
            sigma = this->computeSigma(this->grid->EyPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = -cdt / dx;

            this->splitGrid->eyx[idx] = coeff1 * this->splitGrid->eyx[idx] + coeff2 *
                (this->grid->Bz(i + 1, j, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezx[idx] = coeff1 * this->splitGrid->ezx[idx] - coeff2 *
                (this->grid->By(i + 1, j, k) - this->grid->By(i, j, k));

            // Yee Grid indices: Ex(..., j+1/2, ...), Ez(..., j+1/2, ...)-> sigmaY(j+1/2)
            sigma = this->computeSigma(this->grid->EzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dy);
            else coeff2 = -cdt / dy;

            this->splitGrid->exy[idx] = coeff1 * this->splitGrid->exy[idx] - coeff2 *
                (this->grid->Bz(i, j + 1, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezy[idx] = coeff1 * this->splitGrid->ezy[idx] + coeff2 *
                (this->grid->Bx(i, j + 1, k) - this->grid->Bx(i, j, k));

            // Yee Grid indices: Ex(..., ..., k+1/2), Ey(..., ..., k+1/2)-> sigmaZ(k+1/2)
            sigma = this->computeSigma(this->grid->ExPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dz);
            else coeff2 = -cdt / dz;

            this->splitGrid->exz[idx] = coeff1 * this->splitGrid->exz[idx] + coeff2 *
                (this->grid->By(i, j, k + 1) - this->grid->By(i, j, k));
            this->splitGrid->eyz[idx] = coeff1 * this->splitGrid->eyz[idx] - coeff2 *
                (this->grid->Bx(i, j, k + 1) - this->grid->Bx(i, j, k));

            this->grid->Ex(i, j, k) = this->splitGrid->exy[idx] + this->splitGrid->exz[idx];
            this->grid->Ey(i, j, k) = this->splitGrid->eyx[idx] + this->splitGrid->eyz[idx];
            this->grid->Ez(i, j, k) = this->splitGrid->ezx[idx] + this->splitGrid->ezy[idx];
        }
    }

    inline void PmlFdtd::updateE2D()
    {
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: Ey(i+1/2, ..., ...), Ez(i+1/2, ..., ...) -> sigmaX(i+1/2)
            sigma = this->computeSigma(this->grid->EyPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = -cdt / dx;

            this->splitGrid->eyx[idx] = coeff1 * this->splitGrid->eyx[idx] + coeff2 *
                (this->grid->Bz(i + 1, j, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezx[idx] = coeff1 * this->splitGrid->ezx[idx] - coeff2 *
                (this->grid->By(i + 1, j, k) - this->grid->By(i, j, k));

            // Yee Grid indices: Ex(..., j+1/2, ...), Ez(..., j+1/2, ...)-> sigmaY(j+1/2)
            sigma = this->computeSigma(this->grid->EzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dy);
            else coeff2 = -cdt / dy;

            this->splitGrid->exy[idx] = coeff1 * this->splitGrid->exy[idx] - coeff2 *
                (this->grid->Bz(i, j + 1, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezy[idx] = coeff1 * this->splitGrid->ezy[idx] + coeff2 *
                (this->grid->Bx(i, j + 1, k) - this->grid->Bx(i, j, k));

            // Yee Grid indices: Ex(..., ..., k+1/2), Ey(..., ..., k+1/2)-> sigmaZ(k+1/2)
            sigma = this->computeSigma(this->grid->ExPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);

            this->splitGrid->exz[idx] = coeff1 * this->splitGrid->exz[idx];
            this->splitGrid->eyz[idx] = coeff1 * this->splitGrid->eyz[idx];

            this->grid->Ex(i, j, k) = this->splitGrid->exy[idx] + this->splitGrid->exz[idx];
            this->grid->Ey(i, j, k) = this->splitGrid->eyx[idx] + this->splitGrid->eyz[idx];
            this->grid->Ez(i, j, k) = this->splitGrid->ezx[idx] + this->splitGrid->ezy[idx];
        }
    }

    inline void PmlFdtd::updateE1D()
    {
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;
            FP sigma = 0, coeff1 = 0, coeff2 = 0;

            // Yee Grid indices: Ey(i+1/2, ..., ...), Ez(i+1/2, ..., ...) -> sigmaX(i+1/2)
            sigma = this->computeSigma(this->grid->EyPosition(i, j, k).x, CoordinateEnum::x);

            coeff1 = exp(-sigma * cdt);
            coeff2 = 0;
            if (sigma >= threshold)
                coeff2 = (coeff1 - (FP)1) / (sigma * dx);
            else coeff2 = -cdt / dx;

            this->splitGrid->eyx[idx] = coeff1 * this->splitGrid->eyx[idx] + coeff2 *
                (this->grid->Bz(i + 1, j, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezx[idx] = coeff1 * this->splitGrid->ezx[idx] + coeff2 *
                (this->grid->By(i + 1, j, k) - this->grid->By(i, j, k));

            // Yee Grid indices: Ex(..., j+1/2, ...), Ez(..., j+1/2, ...)-> sigmaY(j+1/2)
            sigma = this->computeSigma(this->grid->EzPosition(i, j, k).y, CoordinateEnum::y);

            coeff1 = exp(-sigma * cdt);

            this->splitGrid->exy[idx] = coeff1 * this->splitGrid->exy[idx];
            this->splitGrid->ezy[idx] = coeff1 * this->splitGrid->ezy[idx];

            // Yee Grid indices: Ex(..., ..., k+1/2), Ey(..., ..., k+1/2)-> sigmaZ(k+1/2)
            sigma = this->computeSigma(this->grid->ExPosition(i, j, k).z, CoordinateEnum::z);

            coeff1 = exp(-sigma * cdt);

            this->splitGrid->exz[idx] = coeff1 * this->splitGrid->exz[idx];
            this->splitGrid->eyz[idx] = coeff1 * this->splitGrid->eyz[idx];

            this->grid->Ex(i, j, k) = this->splitGrid->exy[idx] + this->splitGrid->exz[idx];
            this->grid->Ey(i, j, k) = this->splitGrid->eyx[idx] + this->splitGrid->eyz[idx];
            this->grid->Ez(i, j, k) = this->splitGrid->ezx[idx] + this->splitGrid->ezy[idx];
        }
    }
}
