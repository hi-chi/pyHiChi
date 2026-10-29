#pragma once

#include "Constants.h"
#include "Grid.h"
#include "FormFactor.h"
#include "ParticleArray.h"
#include "Particle.h"
#include <iostream>
#include <vector>

namespace pfc
{
    template<class TGrid, class DerivedClass, class LocalDerivedClass>
    class CurrentDeposition
    {
    public:
        enum class ZeroizeJ {
            NOT_USE_ZEROIZEJ, USE_ZEROIZEJ
        };

        CurrentDeposition(double _dt) : halfDt(0.5 * _dt) {}

        template<class T_Particle>
        void operator()(TGrid* grid, const T_Particle& particle,
            CurrentDeposition::ZeroizeJ UsingZeroizeJ = CurrentDeposition::ZeroizeJ::NOT_USE_ZEROIZEJ) {
            if (UsingZeroizeJ == ZeroizeJ::USE_ZEROIZEJ)
                grid->zeroizeJ();
            static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle);
        }

        template<class T_ParticleArray>
        void operator()(TGrid* grid, T_ParticleArray* particleArray) {
            typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
            grid->zeroizeJ();
            int num_threads = omp_get_max_threads();
            std::vector<LocalDerivedClass> Depositor(num_threads, LocalDerivedClass(Int3(0,0,0),
               static_cast<DerivedClass*>(this), grid));
#pragma omp parallel for
            for (int i = 0; i < particleArray->size(); i++) {
                int thread_num = omp_get_thread_num();
                ParticleProxyType particle = (*particleArray)[i];
                Int3 baseGridIdx = grid->getBaseIndex(
                    particle.getPosition() - (particle.getVelocity() * halfDt));
                if (Depositor[thread_num].getBaseGridIdx() != baseGridIdx) {
                    Depositor[thread_num].addCurrents();
                    Depositor[thread_num].setNewGridCell(baseGridIdx);
                }
                static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle, Depositor[thread_num]);
            }
#pragma omp parallel for
           for (int i = 0; i < num_threads; ++i)
               Depositor[i].addCurrents();
        }

        template<class T_Particle>
        void depositOneParticle(TGrid* grid, T_Particle* particle, LocalDerivedClass& Depositor)
        {
            Depositor.depositCurrent(particle);
        }

        const double halfDt;
    };

    template<class TGrid, class TCurrentDeposition, int BlockSize>
    class LocalDeposition
    {
    public:

        static const int blockSize = BlockSize;

        LocalDeposition(const Int3& _baseGridIdx, TCurrentDeposition* currentDeposition, TGrid* _grid) :
            baseGridIdx(_baseGridIdx), grid(_grid), halfDt(currentDeposition->halfDt)
        {            
            blockOffset = (blockSize) / 2;
            startGridIdx = baseGridIdx - Int3(blockOffset, blockOffset, blockOffset);

            for (int i = 0; i < this->blockSize; ++i)
                for (int j = 0; j < this->blockSize; ++j)
                    for (int k = 0; k < this->blockSize; ++k)
                    {
                        Jx[i][j][k] = 0;
                        Jy[i][j][k] = 0;
                        Jz[i][j][k] = 0;
                    }
        }

        void setNewGridCell(const Int3& _baseGridIdx) {
            baseGridIdx = _baseGridIdx;
            startGridIdx = baseGridIdx - Int3(blockOffset, blockOffset, blockOffset);
        }

        void addCurrents()
        {
            for (int i = 0; i < blockSize; ++i)
                for (int j = 0; j < blockSize; ++j)
                    for (int k = 0; k < blockSize; ++k)
                    {
                        Int3 idx = remainder(startGridIdx + Int3(i, j, k), grid->numCells);
                        #pragma omp atomic
                        grid->Jx(idx) += Jx[i][j][k];
                        #pragma omp atomic
                        grid->Jy(idx) += Jy[i][j][k];
                        #pragma omp atomic
                        grid->Jz(idx) += Jz[i][j][k];
                    }

            for (int i = 0; i < this->blockSize; ++i)
                for (int j = 0; j < this->blockSize; ++j)
                    for (int k = 0; k < this->blockSize; ++k)
                    {
                        Jx[i][j][k] = 0;
                        Jy[i][j][k] = 0;
                        Jz[i][j][k] = 0;
                    }
        }

        template<class TParticle>
        void depositCurrent(const TParticle& particle) {}
        Int3 getBaseGridIdx() { return baseGridIdx; }
    protected:

        FP Jx[blockSize][blockSize][blockSize];
        FP Jy[blockSize][blockSize][blockSize];
        FP Jz[blockSize][blockSize][blockSize];
        TGrid* grid;
        Int3 startGridIdx;
        int blockOffset;
        const double halfDt;
        Int3 baseGridIdx;
    };

    template<class TGrid>
    class CurrentDepositionCIC;

    template<class TGrid>
    class LocalDepositionCIC : public LocalDeposition<TGrid, CurrentDepositionCIC<TGrid>, 3>
    {
    public:
        LocalDepositionCIC(const Int3& blockIdx, CurrentDepositionCIC<TGrid>* currentDeposition, TGrid* _grid)
            : LocalDeposition<TGrid, CurrentDepositionCIC<TGrid>, 3>(blockIdx, currentDeposition, _grid) {}

        template<class T_Particle>
        void depositCurrent(T_Particle* particle)
        {
            FP3 particlePosition = particle->getPosition() - (particle->getVelocity() * this->halfDt);
            FP3 current = (particle->getVelocity() * particle->getCharge() * particle->getWeight()) /
                this->grid->steps.volume();

            Int3 idxJx, idxJy, idxJz;
            FP3 internalCoordsJx, internalCoordsJy, internalCoordsJz;
            idxJx = this->grid->getIndexJx(particlePosition) - this->startGridIdx;
            internalCoordsJx = this->grid->getInternalCoordsJx(particlePosition);

            idxJy = this->grid->getIndexJy(particlePosition) - this->startGridIdx;
            internalCoordsJy = this->grid->getInternalCoordsJy(particlePosition);

            idxJz = this->grid->getIndexJz(particlePosition) - this->startGridIdx;
            internalCoordsJz = this->grid->getInternalCoordsJz(particlePosition);

            FormFactorCIC formFactorJx, formFactorJy, formFactorJz;
            
            formFactorJx(internalCoordsJx);
            formFactorJy(internalCoordsJy);
            formFactorJz(internalCoordsJz);
            depositComponent(idxJx, current.x, formFactorJx, this->Jx);
            depositComponent(idxJy, current.y, formFactorJy, this->Jy);
            depositComponent(idxJz, current.z, formFactorJz, this->Jz);
        }

    private:
        void depositComponent(const Int3 & idx, const FP & value, FormFactorCIC& formFactor,
            FP current[3][3][3])
        {
            for (int i = 0; i <= 1; i++) {
                for (int j = 0; j <= 1; j++) {
                    for (int k = 0; k <= 1; k++) {
                        current[idx.x + i][idx.y + j][idx.z + k] += 
                            (formFactor.c[0][i] *
                             formFactor.c[1][j] *
                             formFactor.c[2][k]
                            ) * value;
                    }
                }
            }
        }
    };

    template<class TGrid>
    class CurrentDepositionCIC : public CurrentDeposition<TGrid, CurrentDepositionCIC<TGrid>,
        LocalDepositionCIC<TGrid>>
    {
    public:
        CurrentDepositionCIC(double _dt) : CurrentDeposition<TGrid, CurrentDepositionCIC<TGrid>,
            LocalDepositionCIC<TGrid>>(_dt) {}
    };

    template<class TGrid>
    class CurrentDepositionTSC;

    template<class TGrid>
    class LocalDepositionTSC : public LocalDeposition<TGrid, CurrentDepositionTSC<TGrid>, 5>
    {
    public:

        LocalDepositionTSC(const Int3 & blockIdx, CurrentDepositionTSC<TGrid>* currentDeposition, TGrid* _grid)
            : LocalDeposition<TGrid, CurrentDepositionTSC<TGrid>, 5>(blockIdx, currentDeposition, _grid) {}

        template<class T_Particle>
        void depositCurrent(T_Particle* particle)
        {
            FP3 particlePosition = particle->getPosition() - (particle->getVelocity() * this->halfDt);
            FP3 current = (particle->getVelocity() * particle->getCharge() * particle->getWeight()) /
                this->grid->steps.volume();

            Int3 idxJx, idxJy, idxJz;
            FP3 internalCoordsJx, internalCoordsJy, internalCoordsJz;
            idxJx = this->grid->getClosestIndexJx(particlePosition) - this->startGridIdx;
            internalCoordsJx = this->grid->getClosestInternalCoordsJx(particlePosition);

            idxJy = this->grid->getClosestIndexJy(particlePosition) - this->startGridIdx;
            internalCoordsJy = this->grid->getClosestInternalCoordsJy(particlePosition);

            idxJz = this->grid->getClosestIndexJz(particlePosition) - this->startGridIdx;
            internalCoordsJz = this->grid->getClosestInternalCoordsJz(particlePosition);

            FormFactorTSC formFactorJx, formFactorJy, formFactorJz;
            
            formFactorJx(internalCoordsJx);
            formFactorJy(internalCoordsJy);
            formFactorJz(internalCoordsJz);
            depositComponent(idxJx, current.x, formFactorJx, this->Jx);
            depositComponent(idxJy, current.y, formFactorJy, this->Jy);
            depositComponent(idxJz, current.z, formFactorJz, this->Jz);
        }

    private:

        void depositComponent(const Int3 & idx,
            const FP & value, FormFactorTSC& formFactor, FP current[5][5][5])
        {
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    for (int k = -1; k <= 1; k++) {
                        current[idx.x + i][idx.y + j][idx.z + k] += 
                            ( formFactor.c[0][i + 1] 
                            * formFactor.c[1][j + 1]
                            * formFactor.c[2][k + 1])
                            * value;
                    }
                }
            }
        }
    };

    template<class TGrid>
    class CurrentDepositionTSC : public CurrentDeposition<TGrid, CurrentDepositionTSC<TGrid>,
        LocalDepositionTSC<TGrid>>
    {
    public:
        CurrentDepositionTSC(double _dt) : CurrentDeposition<TGrid, CurrentDepositionTSC<TGrid>,
            LocalDepositionTSC<TGrid>>(_dt) {}
    };

    template<class TGrid>
    class CurrentDepositionZ1;

    template<class TGrid>
    class LocalDepositionZ1 : public LocalDeposition<TGrid, CurrentDepositionZ1<TGrid>, 3>
    {
    public:

        LocalDepositionZ1(
            const Int3 & blockIdx, CurrentDepositionZ1<TGrid>* currentDeposition, TGrid* _grid)
            : LocalDeposition<TGrid, CurrentDepositionZ1<TGrid>, 3>(blockIdx, currentDeposition, _grid) {}

        template<class T_Particle>
        void depositCurrent(T_Particle* particle)
        {
            FP3 minPosition = this->grid->origin + this->startGridIdx * this->grid->steps;
            FP3 position = particle->getPosition();
            FP3 oldPosition = particle->getPosition() - particle->getVelocity() * (2 * this->halfDt);
            Int3 localIndex = truncate((position - minPosition) / this->grid->steps + FP3(0.5, 0.5, 0.5));
            Int3 oldLocalIndex = truncate((oldPosition - minPosition) / this->grid->steps + FP3(0.5, 0.5, 0.5));
            FP current = particle->getCharge() * particle->getWeight() / this->grid->steps.volume() / (2 * this->halfDt);
            FP3 r;
            for (int i = 0; i < 3; ++i)
            {
                if(oldLocalIndex[i] == localIndex[i]) r[i] = (position[i] + oldPosition[i]) * (FP)0.5;
                else r[i] = (oldLocalIndex[i] + localIndex[i]) * this->grid->steps[i] * (FP)0.5 + minPosition[i];
            }

            FP3 coeff[2], weight1[2], weight2[2];
            weight1[0] = ((oldPosition + r) * (FP)0.5 - minPosition) / this->grid->steps - (FP3)oldLocalIndex + FP3(0.5, 0.5, 0.5);
            weight1[1] = ((position + r) * (FP)0.5 - minPosition) / this->grid->steps - (FP3)localIndex + FP3(0.5, 0.5, 0.5);
            weight2[0] = FP3(1.0, 1.0, 1.0) - weight1[0];
            weight2[1] = FP3(1.0, 1.0, 1.0) - weight1[1];
            coeff[0] = current * (r - oldPosition);
            coeff[1] = current * (position - r);

            Int3 idx[2];
            idx[0] = oldLocalIndex;
            idx[1] = localIndex;

            for(int i = 0; i < 2; i++)
            {
                this->Jx[idx[i].x][idx[i].y]    [idx[i].z]     += coeff[i].x * weight1[i].y * weight1[i].z;
                this->Jx[idx[i].x][idx[i].y - 1][idx[i].z]     += coeff[i].x * weight2[i].y * weight1[i].z;
                this->Jx[idx[i].x][idx[i].y]    [idx[i].z - 1] += coeff[i].x * weight1[i].y * weight2[i].z;
                this->Jx[idx[i].x][idx[i].y - 1][idx[i].z - 1] += coeff[i].x * weight2[i].y * weight2[i].z;

                this->Jy[idx[i].x]    [idx[i].y][idx[i].z]     += coeff[i].y * weight1[i].x * weight1[i].z;
                this->Jy[idx[i].x - 1][idx[i].y][idx[i].z]     += coeff[i].y * weight2[i].x * weight1[i].z;
                this->Jy[idx[i].x]    [idx[i].y][idx[i].z - 1] += coeff[i].y * weight1[i].x * weight2[i].z;
                this->Jy[idx[i].x - 1][idx[i].y][idx[i].z - 1] += coeff[i].y * weight2[i].x * weight2[i].z;

                this->Jz[idx[i].x]    [idx[i].y]    [idx[i].z] += coeff[i].z * weight1[i].x * weight1[i].y;
                this->Jz[idx[i].x - 1][idx[i].y]    [idx[i].z] += coeff[i].z * weight2[i].x * weight1[i].y;
                this->Jz[idx[i].x]    [idx[i].y - 1][idx[i].z] += coeff[i].z * weight1[i].x * weight2[i].y;
                this->Jz[idx[i].x - 1][idx[i].y - 1][idx[i].z] += coeff[i].z * weight2[i].x * weight2[i].y;
            }
        }
    };

    template<class TGrid>
    class CurrentDepositionZ1 : public CurrentDeposition<TGrid, CurrentDepositionZ1<TGrid>,
        LocalDepositionZ1<TGrid>>
    {
    public:
        CurrentDepositionZ1(double _dt) : CurrentDeposition<TGrid, CurrentDepositionZ1<TGrid>,
            LocalDepositionZ1<TGrid>>(_dt) {}
    };

    template<class TGrid>
    class CurrentDepositionZ2;

    template<class TGrid>
    class LocalDepositionZ2 : public LocalDeposition<TGrid, CurrentDepositionZ2<TGrid>, 5>
    {
    public:
        LocalDepositionZ2(
            const Int3 & blockIdx, CurrentDepositionZ2<TGrid>* currentDeposition, TGrid* _grid)
            : LocalDeposition<TGrid, CurrentDepositionZ2<TGrid>, 5>(blockIdx, currentDeposition, _grid) {}

        template<class T_Particle>
        void depositCurrent(T_Particle* particle)
        {
            FP3 position = particle->getPosition();
            FP3 minPosition = this->grid->origin + this->startGridIdx * this->grid->steps;
            
            FP3 oldPosition = position - particle->getVelocity() * (2 * this->halfDt);
            Int3 localIndex = truncate((position - minPosition) / this->grid->steps);
            Int3 oldLocalIndex = truncate((oldPosition - minPosition) / this->grid->steps);
            FP current = particle->getCharge() * particle->getWeight() * particle->getVelocity() * (FP)0.5 
                / this->grid->steps.volume();

            FP3 r;
            for (int i = 0; i < 3; i++)
            {
                if (oldLocalIndex[i] == localIndex[i]) r[i] = (position[i] + oldPosition[i]) * (FP)0.5;
                else r[i] = std::max(oldLocalIndex[i], localIndex[i]) * this->grid->steps[i] + minPosition[i];
            }

            FP3 weight[2];
            weight[0] = ((oldPosition + r) * (FP)0.5 - minPosition) / this->grid->steps - (FP3)oldLocalIndex -
                FP3(0.5, 0.5, 0.5);
            weight[1] = ((position + r) * (FP)0.5 - minPosition) / this->grid->steps - (FP3)localIndex -
                FP3(0.5, 0.5, 0.5);

            FP3 coeff[2][2];
            coeff[0][0] = current * (FP3(0.5, 0.5, 0.5) - weight[0]);
            coeff[1][0] = current * (FP3(0.5, 0.5, 0.5) - weight[1]);
            coeff[0][1] = current * (FP3(0.5, 0.5, 0.5) + weight[0]);
            coeff[1][1] = current * (FP3(0.5, 0.5, 0.5) + weight[1]);

            FP3 weight1[2], weight2[2], weight3[2];
            weight1[0] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) - weight[0]) * (FP3(0.5, 0.5, 0.5) - weight[0]);
            weight1[1] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) - weight[1]) * (FP3(0.5, 0.5, 0.5) - weight[1]);

            weight2[0] = FP3(0.75, 0.75, 0.75) - weight[0] * weight[0];
            weight2[1] = FP3(0.75, 0.75, 0.75) - weight[1] * weight[1];

            weight3[0] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) + weight[0]) * (FP3(0.5, 0.5, 0.5) + weight[0]);
            weight3[1] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) + weight[1]) * (FP3(0.5, 0.5, 0.5) + weight[1]);

            Int3 idx[2];
            idx[0] = oldLocalIndex;
            idx[1] = localIndex;

            for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
            {
                this->Jx[idx[i].x + j][idx[i].y - 1][idx[i].z - 1] += coeff[i][j].x * weight1[i].y * weight1[i].z;
                this->Jx[idx[i].x + j][idx[i].y - 1][idx[i].z]     += coeff[i][j].x * weight1[i].y * weight2[i].z;
                this->Jx[idx[i].x + j][idx[i].y - 1][idx[i].z + 1] += coeff[i][j].x * weight1[i].y * weight3[i].z;
                this->Jx[idx[i].x + j][idx[i].y]    [idx[i].z - 1] += coeff[i][j].x * weight2[i].y * weight1[i].z;
                this->Jx[idx[i].x + j][idx[i].y]    [idx[i].z]     += coeff[i][j].x * weight2[i].y * weight2[i].z;
                this->Jx[idx[i].x + j][idx[i].y]    [idx[i].z + 1] += coeff[i][j].x * weight2[i].y * weight3[i].z;
                this->Jx[idx[i].x + j][idx[i].y + 1][idx[i].z - 1] += coeff[i][j].x * weight3[i].y * weight1[i].z;
                this->Jx[idx[i].x + j][idx[i].y + 1][idx[i].z]     += coeff[i][j].x * weight3[i].y * weight2[i].z;
                this->Jx[idx[i].x + j][idx[i].y + 1][idx[i].z + 1] += coeff[i][j].x * weight3[i].y * weight3[i].z;

                this->Jy[idx[i].x - 1][idx[i].y + j][idx[i].z - 1] += coeff[i][j].y * weight1[i].x * weight1[i].z;
                this->Jy[idx[i].x - 1][idx[i].y + j][idx[i].z]     += coeff[i][j].y * weight1[i].x * weight2[i].z;
                this->Jy[idx[i].x - 1][idx[i].y + j][idx[i].z + 1] += coeff[i][j].y * weight1[i].x * weight3[i].z;
                this->Jy[idx[i].x]    [idx[i].y + j][idx[i].z - 1] += coeff[i][j].y * weight2[i].x * weight1[i].z;
                this->Jy[idx[i].x]    [idx[i].y + j][idx[i].z]     += coeff[i][j].y * weight2[i].x * weight2[i].z;
                this->Jy[idx[i].x]    [idx[i].y + j][idx[i].z + 1] += coeff[i][j].y * weight2[i].x * weight3[i].z;
                this->Jy[idx[i].x + 1][idx[i].y + j][idx[i].z - 1] += coeff[i][j].y * weight3[i].x * weight1[i].z;
                this->Jy[idx[i].x + 1][idx[i].y + j][idx[i].z]     += coeff[i][j].y * weight3[i].x * weight2[i].z;
                this->Jy[idx[i].x + 1][idx[i].y + j][idx[i].z + 1] += coeff[i][j].y * weight3[i].x * weight3[i].z;

                this->Jz[idx[i].x - 1][idx[i].y - 1][idx[i].z + j] += coeff[i][j].z * weight1[i].y * weight1[i].x;
                this->Jz[idx[i].x]    [idx[i].y - 1][idx[i].z + j] += coeff[i][j].z * weight1[i].y * weight2[i].x;
                this->Jz[idx[i].x + 1][idx[i].y - 1][idx[i].z + j] += coeff[i][j].z * weight1[i].y * weight3[i].x;
                this->Jz[idx[i].x - 1][idx[i].y]    [idx[i].z + j] += coeff[i][j].z * weight2[i].y * weight1[i].x;
                this->Jz[idx[i].x]    [idx[i].y]    [idx[i].z + j] += coeff[i][j].z * weight2[i].y * weight2[i].x;
                this->Jz[idx[i].x + 1][idx[i].y]    [idx[i].z + j] += coeff[i][j].z * weight2[i].y * weight3[i].x;
                this->Jz[idx[i].x - 1][idx[i].y + 1][idx[i].z + j] += coeff[i][j].z * weight3[i].y * weight1[i].x;
                this->Jz[idx[i].x]    [idx[i].y + 1][idx[i].z + j] += coeff[i][j].z * weight3[i].y * weight2[i].x;
                this->Jz[idx[i].x + 1][idx[i].y + 1][idx[i].z + j] += coeff[i][j].z * weight3[i].y * weight3[i].x;
            }
        }
    };

    template<class TGrid>
    class CurrentDepositionZ2 : public CurrentDeposition<TGrid, CurrentDepositionZ2<TGrid>,
        LocalDepositionZ2<TGrid>>
    {
    public:
        CurrentDepositionZ2(double _dt) : CurrentDeposition<TGrid, CurrentDepositionZ2<TGrid>,
            LocalDepositionZ2<TGrid>>(_dt) {}
    };

}
