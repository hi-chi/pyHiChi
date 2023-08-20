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
        {
            this->computeCoeffs();
        }

        // constructor for loading
        explicit PmlFdtd(YeeGrid* grid, FP dt, Int3 domainIndexBegin, Int3 domainIndexEnd) :
            PmlReal(grid, dt, domainIndexBegin, domainIndexEnd)
        {}

        void updateB();
        void updateE();

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        // coefficient pre-computing
        void computeCoeffs();

        std::vector<FP> bCoeff1X, bCoeff1Y, bCoeff1Z, eCoeff1X, eCoeff1Y, eCoeff1Z;  // e^(-sigma*dt)
        std::vector<FP> bCoeff2X, bCoeff2Y, bCoeff2Z, eCoeff2X, eCoeff2Y, eCoeff2Z;  // -(e^(-sigma*dt) - 1) / (sigma*dt)

    private:

        void updateB3D();
        void updateB2D();
        void updateB1D();
        void updateE3D();
        void updateE2D();
        void updateE1D();

        void computeCoeffs(
            std::vector<FP>& coeff1X, std::vector<FP>& coeff1Y, std::vector<FP>& coeff1Z,
            std::vector<FP>& coeff2X, std::vector<FP>& coeff2Y, std::vector<FP>& coeff2Z,
            const FP3(YeeGrid::* positionFX)(int, int, int) const,
            const FP3(YeeGrid::* positionFY)(int, int, int) const,
            const FP3(YeeGrid::* positionFZ)(int, int, int) const);
    };

    inline void PmlFdtd::computeCoeffs()
    {
        this->computeCoeffs(bCoeff1X, bCoeff1Y, bCoeff1Z, bCoeff2X, bCoeff2Y, bCoeff2Z,
            &YeeGrid::BxPosition, &YeeGrid::ByPosition, &YeeGrid::BzPosition);
        this->computeCoeffs(eCoeff1X, eCoeff1Y, eCoeff1Z, eCoeff2X, eCoeff2Y, eCoeff2Z,
            &YeeGrid::ExPosition, &YeeGrid::EyPosition, &YeeGrid::EzPosition);
    }

    inline void PmlFdtd::computeCoeffs(
        std::vector<FP>& coeff1X, std::vector<FP>& coeff1Y, std::vector<FP>& coeff1Z,
        std::vector<FP>& coeff2X, std::vector<FP>& coeff2Y, std::vector<FP>& coeff2Z,
        const FP3(YeeGrid::* positionFX)(int, int, int) const,
        const FP3(YeeGrid::* positionFY)(int, int, int) const,
        const FP3(YeeGrid::* positionFZ)(int, int, int) const)
    {
        const int size = this->splitGrid->getNumPmlNodes();

        coeff1X.resize(size);
        coeff1Y.resize(size);
        coeff1Z.resize(size);
        coeff2X.resize(size);
        coeff2Y.resize(size);
        coeff2Z.resize(size);

        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        const FP cdt = constants::c * this->dt;
        const FP threshold = std::numeric_limits<FP>::epsilon();

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index = this->splitGrid->getIndex3d(idx);

            // coordinates according to Yee grid
            FP sigmaX = this->computeSigma((this->grid->*positionFY)(index.x, index.y, index.z).x, CoordinateEnum::x);
            FP sigmaY = this->computeSigma((this->grid->*positionFZ)(index.x, index.y, index.z).y, CoordinateEnum::y);
            FP sigmaZ = this->computeSigma((this->grid->*positionFX)(index.x, index.y, index.z).z, CoordinateEnum::z);

            coeff1X[idx] = exp(-sigmaX * cdt);
            coeff1Y[idx] = exp(-sigmaY * cdt);
            coeff1Z[idx] = exp(-sigmaZ * cdt);

            if (this->grid->dimensionality >= 1)
                coeff2X[idx] = sigmaX >= threshold ? (coeff1X[idx] - (FP)1) / (sigmaX * dx) : -cdt / dx;
            if (this->grid->dimensionality >= 2)
                coeff2Y[idx] = sigmaY >= threshold ? (coeff1Y[idx] - (FP)1) / (sigmaY * dy) : -cdt / dy;
            if (this->grid->dimensionality >= 3)
                coeff2Z[idx] = sigmaZ >= threshold ? (coeff1Z[idx] - (FP)1) / (sigmaZ * dz) : -cdt / dz;
        }
    }

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
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;

            this->splitGrid->byx[idx] = bCoeff1X[idx] * this->splitGrid->byx[idx] - bCoeff2X[idx] *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i - 1, j, k));
            this->splitGrid->bzx[idx] = bCoeff1X[idx] * this->splitGrid->bzx[idx] + bCoeff2X[idx] *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i - 1, j, k));

            this->splitGrid->bxy[idx] = bCoeff1Y[idx] * this->splitGrid->bxy[idx] + bCoeff2Y[idx] *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i, j - 1, k));
            this->splitGrid->bzy[idx] = bCoeff1Y[idx] * this->splitGrid->bzy[idx] - bCoeff2Y[idx] *
                (this->grid->Ex(i, j, k) - this->grid->Ex(i, j - 1, k));

            this->splitGrid->bxz[idx] = bCoeff1Z[idx] * this->splitGrid->bxz[idx] - bCoeff2Z[idx] *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i, j, k - 1));
            this->splitGrid->byz[idx] = bCoeff1Z[idx] * this->splitGrid->byz[idx] + bCoeff2Z[idx] *
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
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;

            this->splitGrid->byx[idx] = bCoeff1X[idx] * this->splitGrid->byx[idx] - bCoeff2X[idx] *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i - 1, j, k));
            this->splitGrid->bzx[idx] = bCoeff1X[idx] * this->splitGrid->bzx[idx] + bCoeff2X[idx] *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i - 1, j, k));

            this->splitGrid->bxy[idx] = bCoeff1Y[idx] * this->splitGrid->bxy[idx] + bCoeff2Y[idx] *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i, j - 1, k));
            this->splitGrid->bzy[idx] = bCoeff1Y[idx] * this->splitGrid->bzy[idx] - bCoeff2Y[idx] *
                (this->grid->Ex(i, j, k) - this->grid->Ex(i, j - 1, k));

            this->splitGrid->bxz[idx] = bCoeff1Z[idx] * this->splitGrid->bxz[idx];
            this->splitGrid->byz[idx] = bCoeff1Z[idx] * this->splitGrid->byz[idx];

            this->grid->Bx(i, j, k) = this->splitGrid->bxy[idx] + this->splitGrid->bxz[idx];
            this->grid->By(i, j, k) = this->splitGrid->byx[idx] + this->splitGrid->byz[idx];
            this->grid->Bz(i, j, k) = this->splitGrid->bzx[idx] + this->splitGrid->bzy[idx];
        }
    }

    inline void PmlFdtd::updateB1D()
    {
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;

            this->splitGrid->byx[idx] = bCoeff1X[idx] * this->splitGrid->byx[idx] - bCoeff2X[idx] *
                (this->grid->Ez(i, j, k) - this->grid->Ez(i - 1, j, k));
            this->splitGrid->bzx[idx] = bCoeff1X[idx] * this->splitGrid->bzx[idx] + bCoeff2X[idx] *
                (this->grid->Ey(i, j, k) - this->grid->Ey(i - 1, j, k));

            this->splitGrid->bxy[idx] = bCoeff1Y[idx] * this->splitGrid->bxy[idx];
            this->splitGrid->bzy[idx] = bCoeff1Y[idx] * this->splitGrid->bzy[idx];

            this->splitGrid->bxz[idx] = bCoeff1Z[idx] * this->splitGrid->bxz[idx];
            this->splitGrid->byz[idx] = bCoeff1Z[idx] * this->splitGrid->byz[idx];

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
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;

            this->splitGrid->eyx[idx] = eCoeff1X[idx] * this->splitGrid->eyx[idx] + eCoeff2X[idx] *
                (this->grid->Bz(i + 1, j, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezx[idx] = eCoeff1X[idx] * this->splitGrid->ezx[idx] - eCoeff2X[idx] *
                (this->grid->By(i + 1, j, k) - this->grid->By(i, j, k));

            this->splitGrid->exy[idx] = eCoeff1Y[idx] * this->splitGrid->exy[idx] - eCoeff2Y[idx] *
                (this->grid->Bz(i, j + 1, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezy[idx] = eCoeff1Y[idx] * this->splitGrid->ezy[idx] + eCoeff2Y[idx] *
                (this->grid->Bx(i, j + 1, k) - this->grid->Bx(i, j, k));

            this->splitGrid->exz[idx] = eCoeff1Z[idx] * this->splitGrid->exz[idx] + eCoeff2Z[idx] *
                (this->grid->By(i, j, k + 1) - this->grid->By(i, j, k));
            this->splitGrid->eyz[idx] = eCoeff1Z[idx] * this->splitGrid->eyz[idx] - eCoeff2Z[idx] *
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
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;

            this->splitGrid->eyx[idx] = eCoeff1X[idx] * this->splitGrid->eyx[idx] + eCoeff2X[idx] *
                (this->grid->Bz(i + 1, j, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezx[idx] = eCoeff1X[idx] * this->splitGrid->ezx[idx] - eCoeff2X[idx] *
                (this->grid->By(i + 1, j, k) - this->grid->By(i, j, k));

            this->splitGrid->exy[idx] = eCoeff1Y[idx] * this->splitGrid->exy[idx] - eCoeff2Y[idx] *
                (this->grid->Bz(i, j + 1, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezy[idx] = eCoeff1Y[idx] * this->splitGrid->ezy[idx] + eCoeff2Y[idx] *
                (this->grid->Bx(i, j + 1, k) - this->grid->Bx(i, j, k));

            this->splitGrid->exz[idx] = eCoeff1Z[idx] * this->splitGrid->exz[idx];
            this->splitGrid->eyz[idx] = eCoeff1Z[idx] * this->splitGrid->eyz[idx];

            this->grid->Ex(i, j, k) = this->splitGrid->exy[idx] + this->splitGrid->exz[idx];
            this->grid->Ey(i, j, k) = this->splitGrid->eyx[idx] + this->splitGrid->eyz[idx];
            this->grid->Ez(i, j, k) = this->splitGrid->ezx[idx] + this->splitGrid->ezy[idx];
        }
    }

    inline void PmlFdtd::updateE1D()
    {
        const int size = this->splitGrid->getNumPmlNodes();
        const FP cdt = constants::c * this->dt;
        const FP dx = this->grid->steps.x, dy = this->grid->steps.y, dz = this->grid->steps.z;

        OMP_FOR()
        for (int idx = 0; idx < size; ++idx)
        {
            Int3 index3d = this->splitGrid->getIndex3d(idx);
            int i = index3d.x, j = index3d.y, k = index3d.z;

            this->splitGrid->eyx[idx] = eCoeff1X[idx] * this->splitGrid->eyx[idx] + eCoeff2X[idx] *
                (this->grid->Bz(i + 1, j, k) - this->grid->Bz(i, j, k));
            this->splitGrid->ezx[idx] = eCoeff1X[idx] * this->splitGrid->ezx[idx] - eCoeff2X[idx] *
                (this->grid->By(i + 1, j, k) - this->grid->By(i, j, k));

            this->splitGrid->exy[idx] = eCoeff1Y[idx] * this->splitGrid->exy[idx];
            this->splitGrid->ezy[idx] = eCoeff1Y[idx] * this->splitGrid->ezy[idx];

            this->splitGrid->exz[idx] = eCoeff1Z[idx] * this->splitGrid->exz[idx];
            this->splitGrid->eyz[idx] = eCoeff1Z[idx] * this->splitGrid->eyz[idx];

            this->grid->Ex(i, j, k) = this->splitGrid->exy[idx] + this->splitGrid->exz[idx];
            this->grid->Ey(i, j, k) = this->splitGrid->eyx[idx] + this->splitGrid->eyz[idx];
            this->grid->Ez(i, j, k) = this->splitGrid->ezx[idx] + this->splitGrid->ezy[idx];
        }
    }

    inline void PmlFdtd::save(std::ostream& ostr)
    {
        PmlReal<YeeGrid>::save(ostr);

        const int size = this->splitGrid->getNumPmlNodes();
        ostr.write((char*)&size, sizeof(size));

        ostr.write((char*)bCoeff1X.data(), sizeof(FP) * size);
        ostr.write((char*)bCoeff1Y.data(), sizeof(FP) * size);
        ostr.write((char*)bCoeff1Z.data(), sizeof(FP) * size);
        ostr.write((char*)bCoeff2X.data(), sizeof(FP) * size);
        ostr.write((char*)bCoeff2Y.data(), sizeof(FP) * size);
        ostr.write((char*)bCoeff2Z.data(), sizeof(FP) * size);

        ostr.write((char*)eCoeff1X.data(), sizeof(FP) * size);
        ostr.write((char*)eCoeff1Y.data(), sizeof(FP) * size);
        ostr.write((char*)eCoeff1Z.data(), sizeof(FP) * size);
        ostr.write((char*)eCoeff2X.data(), sizeof(FP) * size);
        ostr.write((char*)eCoeff2Y.data(), sizeof(FP) * size);
        ostr.write((char*)eCoeff2Z.data(), sizeof(FP) * size);
    }

    inline void PmlFdtd::load(std::istream& istr)
    {
        PmlReal<YeeGrid>::load(istr);

        int size = 0;
        istr.read((char*)&size, sizeof(size));

        bCoeff1X.resize(size); bCoeff1Y.resize(size); bCoeff1Z.resize(size);
        bCoeff2X.resize(size); bCoeff2Y.resize(size); bCoeff2Z.resize(size);

        eCoeff1X.resize(size); eCoeff1Y.resize(size); eCoeff1Z.resize(size);
        eCoeff2X.resize(size); eCoeff2Y.resize(size); eCoeff2Z.resize(size);

        istr.read((char*)bCoeff1X.data(), sizeof(FP) * size);
        istr.read((char*)bCoeff1Y.data(), sizeof(FP) * size);
        istr.read((char*)bCoeff1Z.data(), sizeof(FP) * size);
        istr.read((char*)bCoeff2X.data(), sizeof(FP) * size);
        istr.read((char*)bCoeff2Y.data(), sizeof(FP) * size);
        istr.read((char*)bCoeff2Z.data(), sizeof(FP) * size);

        istr.read((char*)eCoeff1X.data(), sizeof(FP) * size);
        istr.read((char*)eCoeff1Y.data(), sizeof(FP) * size);
        istr.read((char*)eCoeff1Z.data(), sizeof(FP) * size);
        istr.read((char*)eCoeff2X.data(), sizeof(FP) * size);
        istr.read((char*)eCoeff2Y.data(), sizeof(FP) * size);
        istr.read((char*)eCoeff2Z.data(), sizeof(FP) * size);
    }
}
