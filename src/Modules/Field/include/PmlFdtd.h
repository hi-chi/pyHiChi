#pragma once
#include "Grid.h"
#include "FieldSolver.h"
#include "Fdtd.h"
#include "Pml.h"

namespace pfc {
    class FDTD;

    class PmlFdtd : public PmlReal<YeeGridType>
    {
    public:
        PmlFdtd(FDTD * solver, Int3 sizePML) :
            PmlReal((RealFieldSolver<YeeGridType>*)solver, sizePML) {}

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
        if (fieldSolver->grid->dimensionality == 3)
            updateB3D();
        else if (fieldSolver->grid->dimensionality == 2)
            updateB2D();
        else if (fieldSolver->grid->dimensionality == 1)
            updateB1D();
    }


    inline void PmlFdtd::updateB3D()
    {
        // For all cells (i, j, k) in PML use following computational scheme
        // with precomputed coefficients coeffBa, coeffBb:
        //
        // byx(i, j, k) = coeffBa.x(i, j, k) * byx(i, j, k) +
        //     coeffBb.x(i, j, k) * (e.z(i, j, k) - e.z(i - 1, j, k)),
        // bzx(i, j, k) = coeffBa.x(i, j, k) * bzx(i, j, k) -
        //     coeffBb.x(i, j, k) * (e.y(i, j, k) - e.y(i - 1, j, k));
        //
        // bxy(i, j, k) = coeffBa.y(i, j, k) * bxy(i, j, k) -
        //     coeffBb.y(i, j, k) * (e.z(i, j, k) - e.z(i, j - 1, k)),
        // bzy(i, j, k) = coeffBa.y(i, j, k) * bzy(i, j, k) +
        //     coeffBb.y(i, j, k) * (e.x(i, j, k) - e.x(i, j - 1, k));
        //
        // bxz(i, j, k) = coeffBa.z(i, j, k) * bxz(i, j, k) +
        //     coeffBb.z(i, j, k) * (e.y(i, j, k) - e.y(i, j, k - 1)),
        // byz(i, j, k) = coeffBa.z(i, j, k) * byz(i, j, k) -
        //     coeffBb.z(i, j, k) * (e.x(i, j, k) - e.x(i, j, k - 1));
        //
        // b.x(i, j, k) = bxy(i, j, k) + bxz(i, j, k),
        // b.y(i, j, k) = byx(i, j, k) + byz(i, j, k),
        // b.z(i, j, k) = bzx(i, j, k) + bzy(i, j, k).
        YeeGrid * grid = fieldSolver->grid;
#pragma omp parallel for
        for (int idx = 0; idx < numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            byx[idx] = coeffBa[idx].x * byx[idx] + coeffBb[idx].x *
                (grid->Ez(i, j, k) - grid->Ez(i - 1, j, k));
            bzx[idx] = coeffBa[idx].x * bzx[idx] - coeffBb[idx].x *
                (grid->Ey(i, j, k) - grid->Ey(i - 1, j, k));

            bxy[idx] = coeffBa[idx].y * bxy[idx] - coeffBb[idx].y *
                (grid->Ez(i, j, k) - grid->Ez(i, j - 1, k));
            bzy[idx] = coeffBa[idx].y * bzy[idx] + coeffBb[idx].y *
                (grid->Ex(i, j, k) - grid->Ex(i, j - 1, k));

            bxz[idx] = coeffBa[idx].z * bxz[idx] + coeffBb[idx].z *
                (grid->Ey(i, j, k) - grid->Ey(i, j, k - 1));
            byz[idx] = coeffBa[idx].z * byz[idx] - coeffBb[idx].z *
                (grid->Ex(i, j, k) - grid->Ex(i, j, k - 1));

            grid->Bx(i, j, k) = bxy[idx] + bxz[idx];
            grid->By(i, j, k) = byx[idx] + byz[idx];
            grid->Bz(i, j, k) = bzx[idx] + bzy[idx];
        }
    }


    inline void PmlFdtd::updateB2D()
    {
        YeeGrid * grid = fieldSolver->grid;
#pragma omp parallel for
        for (int idx = 0; idx < numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            byx[idx] = coeffBa[idx].x * byx[idx] + coeffBb[idx].x *
                (grid->Ez(i, j, k) - grid->Ez(i - 1, j, k));
            bzx[idx] = coeffBa[idx].x * bzx[idx] - coeffBb[idx].x *
                (grid->Ey(i, j, k) - grid->Ey(i - 1, j, k));

            bxy[idx] = coeffBa[idx].y * bxy[idx] - coeffBb[idx].y *
                (grid->Ez(i, j, k) - grid->Ez(i, j - 1, k));
            bzy[idx] = coeffBa[idx].y * bzy[idx] + coeffBb[idx].y *
                (grid->Ex(i, j, k) - grid->Ex(i, j - 1, k));

            bxz[idx] = coeffBa[idx].z * bxz[idx];
            byz[idx] = coeffBa[idx].z * byz[idx];

            grid->Bx(i, j, k) = bxy[idx] + bxz[idx];
            grid->By(i, j, k) = byx[idx] + byz[idx];
            grid->Bz(i, j, k) = bzx[idx] + bzy[idx];
        }
    }


    inline void PmlFdtd::updateB1D()
    {
        YeeGrid * grid = fieldSolver->grid;
#pragma omp parallel for
        for (int idx = 0; idx < numCells; ++idx)
        {
            int i = cellIndex[idx].x;
            int j = cellIndex[idx].y;
            int k = cellIndex[idx].z;

            byx[idx] = coeffBa[idx].x * byx[idx] + coeffBb[idx].x *
                (grid->Ez(i, j, k) - grid->Ez(i - 1, j, k));
            bzx[idx] = coeffBa[idx].x * bzx[idx] - coeffBb[idx].x *
                (grid->Ey(i, j, k) - grid->Ey(i - 1, j, k));

            bxy[idx] = coeffBa[idx].y * bxy[idx];
            bzy[idx] = coeffBa[idx].y * bzy[idx];

            bxz[idx] = coeffBa[idx].z * bxz[idx];
            byz[idx] = coeffBa[idx].z * byz[idx];

            grid->Bx(i, j, k) = bxy[idx] + bxz[idx];
            grid->By(i, j, k) = byx[idx] + byz[idx];
            grid->Bz(i, j, k) = bzx[idx] + bzy[idx];
        }
    }


    inline void PmlFdtd::updateE()
    {
        if (fieldSolver->grid->dimensionality == 3)
            updateE3D();
        else if (fieldSolver->grid->dimensionality == 2)
            updateE2D();
        else if (fieldSolver->grid->dimensionality == 1)
            updateE1D();
    }


    inline void PmlFdtd::updateE3D()
    {
        // For all nodes (i, j, k) in PML use following computational scheme
        // with precomputed coefficients coeffEa, coeffEb:
        //
        // eyx(i, j, k) = coeffEa.x(i, j, k) * eyx(i, j, k) +
        //     coeffEb.x(i, j, k) * (b.z(i + 1, j, k) - b.z(i, j, k)),
        // ezx(i, j, k) = coeffEa.x(i, j, k) * ezx(i, j, k) -
        //     coeffEb.x(i, j, k) * (b.y(i + 1, j, k) - b.y(i, j, k));
        //
        // exy(i, j, k) = coeffEa.y(i, j, k) * exy(i, j, k) -
        //     coeffEb.y(i, j, k) * (b.z(i, j + 1, k) - b.z(i, j, k)),
        // ezy(i, j, k) = coeffEa.y(i, j, k) * ezy(i, j, k) +
        //     coeffEb.y(i, j, k) * (b.x(i, j + 1, k) - b.x(i, j, k));
        //
        // exz(i, j, k) = coeffEa.z(i, j, k) * exz(i, j, k) +
        //     coeffEb.z(i, j, k) * (b.y(i, j, k + 1) - b.y(i, j, k)),
        // eyz(i, j, k) = coeffEa.z(i, j, k) * eyz(i, j, k) -
        //     coeffEb.z(i, j, k) * (b.x(i, j, k + 1) - b.x(i, j, k));
        //
        // e.x(i, j, k) = exy(i, j, k) + exz(i, j, k),
        // e.y(i, j, k) = eyx(i, j, k) + eyz(i, j, k),
        // e.z(i, j, k) = ezx(i, j, k) + ezy(i, j, k).
        YeeGrid * grid = fieldSolver->grid;
        Int3 edgeIdx = grid->numCells - Int3(1, 1, 1);
#pragma omp parallel for
        for (int idx = 0; idx < numNodes; ++idx)
        {
            int i = nodeIndex[idx].x;
            int j = nodeIndex[idx].y;
            int k = nodeIndex[idx].z;

            if (i != edgeIdx.x)
            {
                eyx[idx] = coeffJ[idx] * grid->Jy(i, j, k) + coeffEa[idx].x * eyx[idx] +
                    coeffEb[idx].x * (grid->Bz(i + 1, j, k) - grid->Bz(i, j, k));
                ezx[idx] = coeffJ[idx] * grid->Jz(i, j, k) + coeffEa[idx].x * ezx[idx] -
                    coeffEb[idx].x * (grid->By(i + 1, j, k) - grid->By(i, j, k));
            }

            if (j != edgeIdx.y)
            {
                exy[idx] = coeffJ[idx] * grid->Jx(i, j, k) + coeffEa[idx].y * exy[idx] -
                    coeffEb[idx].y * (grid->Bz(i, j + 1, k) - grid->Bz(i, j, k));
                ezy[idx] = coeffJ[idx] * grid->Jz(i, j, k) + coeffEa[idx].y * ezy[idx] +
                    coeffEb[idx].y * (grid->Bx(i, j + 1, k) - grid->Bx(i, j, k));
            }

            if (k != edgeIdx.z)
            {
                exz[idx] = coeffJ[idx] * grid->Jx(i, j, k) + coeffEa[idx].z * exz[idx] +
                    coeffEb[idx].z * (grid->By(i, j, k + 1) - grid->By(i, j, k));
                eyz[idx] = coeffJ[idx] * grid->Jy(i, j, k) + coeffEa[idx].z * eyz[idx] -
                    coeffEb[idx].z * (grid->Bx(i, j, k + 1) - grid->Bx(i, j, k));
            }

            grid->Ex(i, j, k) = exy[idx] + exz[idx];
            grid->Ey(i, j, k) = eyx[idx] + eyz[idx];
            grid->Ez(i, j, k) = ezx[idx] + ezy[idx];
        }
    }


    inline void PmlFdtd::updateE2D()
    {
        YeeGrid * grid = fieldSolver->grid;
        Int3 edgeIdx = grid->numCells - Int3(1, 1, 1);
#pragma omp parallel for
        for (int idx = 0; idx < numNodes; ++idx)
        {
            int i = nodeIndex[idx].x;
            int j = nodeIndex[idx].y;
            int k = nodeIndex[idx].z;

            if (i != edgeIdx.x)
            {
                eyx[idx] = coeffJ[idx] * grid->Jy(i, j, k) + coeffEa[idx].x * eyx[idx] +
                    coeffEb[idx].x * (grid->Bz(i + 1, j, k) - grid->Bz(i, j, k));
                ezx[idx] = coeffJ[idx] * grid->Jz(i, j, k) + coeffEa[idx].x * ezx[idx] -
                    coeffEb[idx].x * (grid->By(i + 1, j, k) - grid->By(i, j, k));
            }

            if (j != edgeIdx.y)
            {
                exy[idx] = coeffJ[idx] * grid->Jx(i, j, k) + coeffEa[idx].y * exy[idx] -
                    coeffEb[idx].y * (grid->Bz(i, j + 1, k) - grid->Bz(i, j, k));
                ezy[idx] = coeffJ[idx] * grid->Jz(i, j, k) + coeffEa[idx].y * ezy[idx] +
                    coeffEb[idx].y * (grid->Bx(i, j + 1, k) - grid->Bx(i, j, k));
            }

            exz[idx] = coeffJ[idx] * grid->Jx(i, j, k) + coeffEa[idx].z * exz[idx];
            eyz[idx] = coeffJ[idx] * grid->Jy(i, j, k) + coeffEa[idx].z * eyz[idx];

            grid->Ex(i, j, k) = exy[idx] + exz[idx];
            grid->Ey(i, j, k) = eyx[idx] + eyz[idx];
            grid->Ez(i, j, k) = ezx[idx] + ezy[idx];
        }
    }


    inline void PmlFdtd::updateE1D()
    {
        YeeGrid * grid = fieldSolver->grid;
        Int3 edgeIdx = grid->numCells - Int3(1, 1, 1);
#pragma omp parallel for
        for (int idx = 0; idx < numNodes; ++idx)
        {
            int i = nodeIndex[idx].x;
            int j = nodeIndex[idx].y;
            int k = nodeIndex[idx].z;

            if (i != edgeIdx.x)
            {
                eyx[idx] = coeffJ[idx] * grid->Jy(i, j, k) + coeffEa[idx].x * eyx[idx] +
                    coeffEb[idx].x * (grid->Bz(i + 1, j, k) - grid->Bz(i, j, k));
                ezx[idx] = coeffJ[idx] * grid->Jz(i, j, k) + coeffEa[idx].x * ezx[idx] -
                    coeffEb[idx].x * (grid->By(i + 1, j, k) - grid->By(i, j, k));
            }

            exy[idx] = coeffJ[idx] * grid->Jx(i, j, k) + coeffEa[idx].y * exy[idx];
            ezy[idx] = coeffJ[idx] * grid->Jz(i, j, k) + coeffEa[idx].y * ezy[idx];

            exz[idx] = coeffJ[idx] * grid->Jx(i, j, k) + coeffEa[idx].z * exz[idx];
            eyz[idx] = coeffJ[idx] * grid->Jy(i, j, k) + coeffEa[idx].z * eyz[idx];

            grid->Ex(i, j, k) = exy[idx] + exz[idx];
            grid->Ey(i, j, k) = eyx[idx] + eyz[idx];
            grid->Ez(i, j, k) = ezx[idx] + ezy[idx];
        }
    }
}