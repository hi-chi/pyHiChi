#pragma once

#include <random>
#include "Dimension.h"
#include "Particle.h"
#include "VectorsProxy.h"

#include <fstream>


namespace pfc {
    class ParticleGenerator {
        typedef typename VectorTypeHelper<Three, Real>::Type PositionType;
        typedef typename VectorTypeHelper<Three, Real>::Type MomentumType;
        typedef ParticleTypes TypeIndexType;
        typedef Real WeightType;
    public:

        ParticleGenerator():
            genRound(0), genPosition(0), genMomentum(0) {
            // TODO: seed initialization according to the current domain
        }

        template<class T_ParticleArray, GridTypes gridType>
        void operator()(T_ParticleArray* particleArray,
            const Grid<FP, gridType>* grid,
            FP(*particleDensity)(FP, FP, FP),
            FP(*initialTemperature)(FP, FP, FP),
            FP3(*initialMomentum)(FP, FP, FP),
            WeightType weight = 1.0,
            TypeIndexType typeIndex = ParticleTypes::Electron)
        {
            Int3 startIndex = grid->getNumExternalLeftCells();
            Int3 endIndex = grid->getNumExternalLeftCells() + grid->numInternalCells;

            for (int i = startIndex.x; i < endIndex.x; ++i)
                for (int j = startIndex.y; j < endIndex.y; ++j)
                    for (int k = startIndex.z; k < endIndex.z; ++k) {
                        Int3 minIdx(i, j, k); Int3 maxIdx(i + 1, j + 1, k + 1);
                        FP3 minCellCoords = grid->origin + minIdx * grid->steps;
                        FP3 maxCellCoords = grid->origin + maxIdx * grid->steps;

                        this->operator()(particleArray,
                            minCellCoords, maxCellCoords,
                            particleDensity,
                            initialTemperature,
                            initialMomentum,
                            weight
                            );
                    }
        }

        template<class T_ParticleArray>
        void operator()(T_ParticleArray* particleArray,
            const FP3& minCellCoords, const FP3& maxCellCoords,
            FP(*particleDensity)(FP, FP, FP),
            FP(*initialTemperature)(FP, FP, FP),
            FP3(*initialMomentum)(FP, FP, FP),
            WeightType weight = 1.0,
            TypeIndexType typeIndex = ParticleTypes::Electron)
        {
            FP3 center = (minCellCoords + maxCellCoords) * 0.5;
            // number of particle in cell according to density
            FP expectedParticleNum = particleDensity(center.x, center.y, center.z) *
                (maxCellCoords - minCellCoords).volume() / weight;

            int particleNum = int(expectedParticleNum);
            // random shift of expectedParticleNum by 1 to eliminate the effect of rounding
            std::uniform_real_distribution<FP> distRound(0.0, 1.0);
            if (distRound(genRound) < expectedParticleNum - (FP)particleNum)
                ++particleNum;
            for (int i = 0; i < particleNum; ++i) {
                PositionType particlePosition(getParticleRandomPosition(minCellCoords, maxCellCoords));

                MomentumType meanMomentum = initialMomentum(particlePosition.x, particlePosition.y, particlePosition.z);
                FP temperature = initialTemperature(particlePosition.x, particlePosition.y, particlePosition.z);
                MomentumType particleMomentum(getParticleRandomMomentum(meanMomentum, temperature));

                Particle3d newParticle(particlePosition, particleMomentum, weight, typeIndex);
                particleArray->pushBack(newParticle);
            }
        }

    private:
        
        default_random_engine genRound, genPosition, genMomentum;

        FP3 getParticleRandomPosition(const FP3& minCoord, const FP3& maxCoord) {
            std::uniform_real_distribution<FP> distPosX(minCoord.x, maxCoord.x);
            std::uniform_real_distribution<FP> distPosY(minCoord.y, maxCoord.y);
            std::uniform_real_distribution<FP> distPosZ(minCoord.z, maxCoord.z);
            FP3 position(distPosX(genPosition), distPosY(genPosition), distPosZ(genPosition));
            return position;
        }

        MomentumType getParticleRandomMomentum(const MomentumType& meanMomentum, FP temperature) {
            FP alpha = temperature / (constants::electronMass * constants::c * constants::c) + 1.0;
            alpha = 1.5 / (alpha * alpha - 1.0);
            FP sigma = std::sqrt(0.5 / alpha) * constants::electronMass * constants::c;
            std::normal_distribution<FP> distMom(0.0, 1.0);
            MomentumType momentum(
                meanMomentum.x + sigma * distMom(genMomentum),
                meanMomentum.y + sigma * distMom(genMomentum),
                meanMomentum.z + sigma * distMom(genMomentum)
            );
            return momentum;
        }
    };
}
